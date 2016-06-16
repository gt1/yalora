/*
    yalora
    Copyright (C) 2016 German Tischler

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#if 0
#define LIBMAUS2_CHAIN_LINK_DEBUG
#define CHAINNODEINFOSET_DOT
#endif

#include <config.h>

#include <libmaus2/lcs/SMEMProcessor.hpp>

#include <libmaus2/fastx/DNAIndexMetaDataBigBandBiDir.hpp>
#include <libmaus2/fastx/FastAToCompact4BigBandBiDir.hpp>
#include <libmaus2/util/ArgInfo.hpp>
#include <libmaus2/util/ArgParser.hpp>
#include <libmaus2/suffixsort/bwtb3m/BwtMergeSortOptions.hpp>
#include <libmaus2/suffixsort/bwtb3m/BwtMergeSort.hpp>
#include <libmaus2/fastx/StreamFastAReader.hpp>
#include <libmaus2/bambam/BamDecoder.hpp>
#include <libmaus2/bambam/CigarStringParser.hpp>
#include <libmaus2/bambam/BamAlignmentEncoderBase.hpp>
#include <libmaus2/bambam/BamWriter.hpp>
#include <libmaus2/fastx/FastAIndex.hpp>
#include <libmaus2/parallel/SynchronousCounter.hpp>
#include <libmaus2/parallel/LockedGrowingFreeList.hpp>
#include <libmaus2/fastx/FastAIndexGenerator.hpp>
#include <libmaus2/fastx/FastAStreamSet.hpp>
#include <libmaus2/bambam/BamBlockWriterBaseFactory.hpp>
#include <libmaus2/util/MemoryStatistics.hpp>

typedef ::libmaus2::fastx::EntityBuffer<uint8_t,libmaus2::bambam::BamAlignment::D_array_alloc_type> buffer_type;

struct BufferTypeInfo
{
	typedef buffer_type element_type;
	typedef element_type::shared_ptr_type pointer_type;

	static pointer_type getNullPointer()
	{
		return pointer_type();
	}

	static pointer_type deallocate(pointer_type /* p */)
	{
		return getNullPointer();
	}
};

struct BufferAllocator
{
	typedef buffer_type element_type;
	typedef element_type::shared_ptr_type pointer_type;

	pointer_type operator()() const
	{
		return pointer_type(new element_type);
	}
};

struct LockedWriterHeapElement
{
	uint64_t id;
	uint64_t subid;
	buffer_type::shared_ptr_type buffer;

	LockedWriterHeapElement() {}
	LockedWriterHeapElement(std::istream & in, buffer_type::shared_ptr_type rbuffer) : buffer(rbuffer)
	{
		deserialise(in);
	}
	LockedWriterHeapElement(uint64_t const rid, uint64_t const rsubid, buffer_type::shared_ptr_type rbuffer)
	: id(rid), subid(rsubid), buffer(rbuffer) {}

	bool operator<(LockedWriterHeapElement const & other) const
	{
		if ( id != other.id )
			return id < other.id;
		else
			return subid < other.subid;
	}

	void serialise(std::ostream & out) const
	{
		libmaus2::util::NumberSerialisation::serialiseNumber(out,id);
		libmaus2::util::NumberSerialisation::serialiseNumber(out,subid);
		buffer->serialise(out);
	}

	void deserialise(std::istream & in)
	{
		id = libmaus2::util::NumberSerialisation::deserialiseNumber(in);
		subid = libmaus2::util::NumberSerialisation::deserialiseNumber(in);
		buffer->deserialise(in);
	}
};

struct LockedWriterOverflowFile
{
	typedef LockedWriterOverflowFile this_type;
	typedef libmaus2::util::unique_ptr<this_type>::type unique_ptr_type;
	typedef libmaus2::util::shared_ptr<this_type>::type shared_ptr_type;

	std::string const fn;
	libmaus2::aio::OutputStreamInstance::unique_ptr_type OSI;
	uint64_t n;

	LockedWriterOverflowFile(std::string const & rfn)
	: fn(rfn), OSI(new libmaus2::aio::OutputStreamInstance(fn)), n(0)
	{

	}

	void put(LockedWriterHeapElement const & E, libmaus2::parallel::LockedGrowingFreeList<buffer_type,BufferAllocator,BufferTypeInfo> & bufferFreeList)
	{
		E.serialise(*OSI);
		n += 1;
		bufferFreeList.put(E.buffer);
	}

	void flush(
		libmaus2::parallel::LockedGrowingFreeList<buffer_type,BufferAllocator,BufferTypeInfo> & bufferFreeList,
		libmaus2::bambam::BamBlockWriterBase & writer
	)
	{
		OSI->flush();
		OSI.reset();

		libmaus2::aio::InputStreamInstance::unique_ptr_type ISI(new libmaus2::aio::InputStreamInstance(fn));
		libmaus2::autoarray::AutoArray < LockedWriterHeapElement > A(n);
		for ( uint64_t i = 0; i < n; ++i )
			A[i] = LockedWriterHeapElement(*ISI,bufferFreeList.get());
		ISI.reset();

		libmaus2::aio::FileRemoval::removeFile(fn);

		std::sort(A.begin(),A.end());

		for ( uint64_t i = 0; i < n; ++i )
		{
			LockedWriterHeapElement LWHE = A[i];

			buffer_type::shared_ptr_type buffer = LWHE.buffer;

			writer.writeBamBlock(buffer->buffer,buffer->length);

			bufferFreeList.put(buffer);
		}
	}
};

struct LockedWriter
{
	libmaus2::parallel::LockedGrowingFreeList<buffer_type,BufferAllocator,BufferTypeInfo> bufferFreeList;

	libmaus2::bambam::BamBlockWriterBase & writer;
	libmaus2::parallel::PosixSpinLock writerlock;

	uint64_t volatile nextfinish;
	std::set < uint64_t > finishpending;
	libmaus2::parallel::PosixSpinLock finishlock;

	libmaus2::util::FiniteSizeHeap<LockedWriterHeapElement> H;
	libmaus2::parallel::PosixSpinLock hlock;

	uint64_t const overflowsize;

	std::map < uint64_t, LockedWriterOverflowFile::shared_ptr_type > Mexp;
	libmaus2::parallel::PosixSpinLock exlock;

	std::string const tmpprefix;

	libmaus2::util::SimpleQueue < LockedWriterHeapElement > SQ;

	uint64_t getExpungeId(uint64_t const id) const
	{
		return id / overflowsize;
	}

	LockedWriter(
		libmaus2::bambam::BamBlockWriterBase & rwriter,
		std::string const & rtmpprefix)
	: writer(rwriter), nextfinish(0), H(16*1024), overflowsize(128), tmpprefix(rtmpprefix) {}

	void put(uint8_t const * u, uint64_t const l)
	{
		{
			libmaus2::parallel::ScopePosixSpinLock slock(writerlock);
			writer.writeBamBlock(u,l);
		}
	}

	void putex(LockedWriterHeapElement & LWHE)
	{
		uint64_t const exid = getExpungeId(LWHE.id);

		{
			libmaus2::parallel::ScopePosixSpinLock slock(exlock);
			if ( Mexp.find(exid) == Mexp.end() )
			{
				std::ostringstream fnostr;
				fnostr << tmpprefix << "_" << exid;
				std::string const fn = fnostr.str();

				libmaus2::util::TempFileRemovalContainer::addTempFile(fn);

				Mexp[exid] = LockedWriterOverflowFile::shared_ptr_type(new LockedWriterOverflowFile(fn));
			}
			Mexp.find(exid)->second->put(LWHE,bufferFreeList);
		}
	}

	bool isExpunged(uint64_t const id)
	{
		bool r = false;
		uint64_t const exid = getExpungeId(id);

		{
			libmaus2::parallel::ScopePosixSpinLock slock(exlock);
			r = Mexp.find(exid) != Mexp.end();
		}

		return r;
	}

	bool isExpunged(LockedWriterHeapElement const & LWHE)
	{
		return isExpunged(LWHE.id);
	}

	bool getFinishedPending(uint64_t & pending)
	{
		libmaus2::parallel::ScopePosixSpinLock slock(finishlock);
		if ( finishpending.size() && *(finishpending.begin()) == nextfinish )
		{
			pending = nextfinish;
			finishpending.erase(pending);
			return true;
		}
		else
			return false;
	}

	void bumpFinishPending()
	{
		libmaus2::parallel::ScopePosixSpinLock slock(finishlock);
		nextfinish++;
	}

	void flushExpungeFile(LockedWriterOverflowFile::shared_ptr_type Pfile)
	{
		libmaus2::parallel::ScopePosixSpinLock slock(writerlock);
		Pfile->flush(bufferFreeList,writer);
	}

	void put(libmaus2::util::SimpleQueue<buffer_type::shared_ptr_type> & SQB, uint64_t const rid)
	{
		{
			// lock heap
			libmaus2::parallel::ScopePosixSpinLock slock(hlock);

			// if block is already in expunge list
			if ( isExpunged(rid) )
			{
				for ( uint64_t subid = 0; ! SQB.empty(); ++subid )
				{
					buffer_type::shared_ptr_type pbuffer = SQB.pop_front();
					LockedWriterHeapElement LWHE(rid,subid,pbuffer);
					putex(LWHE);
				}
			}
			else
			{
				uint64_t const toinsert = SQB.size();

				while ( (!H.empty()) && toinsert > H.free() )
				{
					uint64_t const minex = getExpungeId(H.top().id);

					while ( !H.empty() && getExpungeId(H.top().id) == minex )
					{
						LockedWriterHeapElement LWHE = H.pop();
						putex(LWHE);
					}
				}

				assert ( H.empty() || toinsert <= H.free() );

				// check whether it fits in heap
				if (
					(!isExpunged(rid))
					&&
					toinsert <= H.free()
				)
				{
					for ( uint64_t subid = 0; ! SQB.empty(); ++subid )
					{
						buffer_type::shared_ptr_type pbuffer = SQB.pop_front();
						LockedWriterHeapElement LWHE(rid,subid,pbuffer);
						assert ( ! H.full() );
						H.push(LWHE);
					}
				}
				// no, expunge it
				else
				{
					assert (
						H.empty()
						||
						isExpunged(rid)
					);

					for ( uint64_t subid = 0; ! SQB.empty(); ++subid )
					{
						buffer_type::shared_ptr_type pbuffer = SQB.pop_front();
						LockedWriterHeapElement LWHE(rid,subid,pbuffer);
						putex(LWHE);
					}
				}
			}

			assert ( SQB.size() == 0 );
		}

		{
			libmaus2::parallel::ScopePosixSpinLock slock(finishlock);
			finishpending.insert(rid);
		}

		uint64_t pending;
		while ( getFinishedPending(pending) )
		{
			LockedWriterOverflowFile::shared_ptr_type Pfile;

			{
				libmaus2::parallel::ScopePosixSpinLock shlock(hlock);

				// is this the last id in the block?
				if ( (pending + 1) % overflowsize == 0 )
				{
					uint64_t const exid = getExpungeId(pending);
					assert ( H.empty() || getExpungeId(H.top().id) >= exid );

					{
						libmaus2::parallel::ScopePosixSpinLock slock(exlock);
						assert ( Mexp.empty() || Mexp.begin()->first >= exid );
					}

					if ( isExpunged(pending) )
					{
						bool const ok = H.empty() || getExpungeId(H.top().id) > exid;
						if ( ! ok )
						{
							assert ( ! H.empty() );

							std::cerr << "exid=" << exid << " top " << H.top().id << " " << getExpungeId(H.top().id) << std::endl;

							assert ( ok );
						}

						libmaus2::parallel::ScopePosixSpinLock slock(exlock);
						std::map < uint64_t, LockedWriterOverflowFile::shared_ptr_type >::iterator it = Mexp.find(exid);

						assert ( it != Mexp.end() );

						Pfile = it->second;
						Mexp.erase(it);
					}
					else
					{
						while ( (!H.empty()) && getExpungeId(H.top().id) == exid )
							SQ.push_back(H.pop());
					}
				}
			}

			if ( (pending + 1) % overflowsize == 0 )
			{
				while ( ! SQ.empty() )
				{
					LockedWriterHeapElement LWHE = SQ.pop_front();
					buffer_type::shared_ptr_type pbuffer = LWHE.buffer;
					put(pbuffer->buffer,pbuffer->length);
					bufferFreeList.put(pbuffer);
				}

				if ( Pfile )
					flushExpungeFile(Pfile);
			}

			bumpFinishPending();
		}
	}

	void flush()
	{
		assert ( Mexp.size() <= 1 );

		while ( !Mexp.empty() )
		{
			std::map < uint64_t, LockedWriterOverflowFile::shared_ptr_type >::iterator it = Mexp.begin();
			LockedWriterOverflowFile::shared_ptr_type Pfile = it->second;
			flushExpungeFile(Pfile);
			Mexp.erase(it);
		}

		while ( !(H.empty()) )
		{
			LockedWriterHeapElement LWHE = H.pop();
			buffer_type::shared_ptr_type pbuffer = LWHE.buffer;
			put(pbuffer->buffer,pbuffer->length);
			bufferFreeList.put(pbuffer);
		}
	}
};

struct LockedGetInterface
{
	typedef LockedGetInterface this_type;
	typedef libmaus2::util::unique_ptr<this_type>::type unique_ptr_type;

	virtual ~LockedGetInterface() {}

	virtual bool getNext(
		std::string & sid,
		std::string & rpat,
		uint64_t & rid,
		libmaus2::bambam::BamAlignment & ralgn,
		libmaus2::bambam::BamAlignment * & palgn
	) = 0;
};

struct BamLockedGet : public LockedGetInterface
{
	libmaus2::parallel::PosixSpinLock lock;
	libmaus2::bambam::BamDecoder dec;
	libmaus2::bambam::BamAlignment & algn;
	uint64_t volatile id;

	BamLockedGet(std::istream & rin) : dec(rin), algn(dec.getAlignment()), id(0)
	{
	}

	bool getNext(
		std::string & sid,
		std::string & rpat,
		uint64_t & rid,
		libmaus2::bambam::BamAlignment & ralgn,
		libmaus2::bambam::BamAlignment * & palgn
	)
	{
		bool ok = false;
		uint64_t lid = 0;
		{
			libmaus2::parallel::ScopePosixSpinLock slock(lock);
			ok = dec.readAlignment();

			if ( ok )
			{
				lid = id++;
				ralgn.swap(algn);
			}
		}

		if ( ok )
		{
			rid = lid;
			sid = ralgn.getName();
			rpat = ralgn.isReverse() ? ralgn.getReadRC() : ralgn.getRead();
			palgn = &ralgn;
		}

		return ok;
	}
};

struct FastALockedGet : public LockedGetInterface
{
	libmaus2::parallel::PosixSpinLock lock;
	libmaus2::fastx::StreamFastAReaderWrapper SFARW;
	libmaus2::fastx::StreamFastAReaderWrapper::pattern_type pattern;
	uint64_t volatile id;

	FastALockedGet(std::istream & rin) : SFARW(rin), id(0)
	{
	}

	bool getNext(
		std::string & sid,
		std::string & rpat,
		uint64_t & rid,
		libmaus2::bambam::BamAlignment & /* ralgn */,
		libmaus2::bambam::BamAlignment * & palgn
	)
	{
		bool ok = false;
		uint64_t lid = 0;
		{
			libmaus2::parallel::ScopePosixSpinLock slock(lock);
			ok = SFARW.getNextPatternUnlocked(pattern);
			if ( ok )
			{
				lid = id++;
				sid = pattern.getShortStringId();
				rpat = pattern.spattern;
			}
		}
		rid = lid;

		palgn = 0;

		return ok;
	}
};

template<typename _ssa_type>
struct AlignContext
{
	typedef _ssa_type ssa_type;
	typedef AlignContext<ssa_type> this_type;
	typedef typename libmaus2::util::unique_ptr<this_type>::type unique_ptr_type;

	LockedWriter & LW;
	libmaus2::fastx::FastAIndex PFAI;
	libmaus2::fastx::DNAIndexMetaDataBigBandBiDir meta;
	libmaus2::rank::DNARank const & Prank;
	libmaus2::rank::DNARankKmerCache const & Pcache;
	libmaus2::fastx::CoordinateCacheBiDir cocache;
	uint64_t n;
	char const * text;
	uint64_t const minfreq;
	uint64_t const minlen;
	uint64_t const limit;
	uint64_t const minsplitlength;
	uint64_t const minsplitsize;
	libmaus2::lcs::SMEMProcessor<ssa_type> proc;
	libmaus2::lcs::NNP nnp;
	libmaus2::lcs::NNPTraceContainer trace;
	libmaus2::lcs::AlignmentTraceContainer ATC;
	::libmaus2::fastx::EntityBuffer<uint8_t,libmaus2::bambam::BamAlignment::D_array_alloc_type> buffer;
	uint64_t onrefid;
	uint64_t on;
	uint64_t off;
	libmaus2::parallel::SynchronousCounter<uint64_t> & cnt;
	libmaus2::timing::RealTimeClock & rtc;
	libmaus2::autoarray::AutoArray<char> Hbuf;
	libmaus2::util::SimpleQueue < buffer_type::shared_ptr_type > SQB;

	libmaus2::autoarray::AutoArray< std::pair<libmaus2::lcs::AlignmentTraceContainer::step_type,uint64_t> > Aopblocks;
	libmaus2::autoarray::AutoArray< libmaus2::bambam::cigar_operation> Aop;

	uint64_t minalgnlen;

	AlignContext(
		LockedWriter & rLW,
		libmaus2::fastx::FastAIndex const & rPFAI,
		libmaus2::fastx::DNAIndexMetaDataBigBandBiDir & rmeta,
		libmaus2::rank::DNARank const & rPrank,
		libmaus2::rank::DNARankKmerCache const & rPcache,
		ssa_type const & rBSSSA,
		char const * rtext,
		uint64_t const rminfreq,
		uint64_t const rminlen,
		uint64_t const rlimit,
		uint64_t const rminsplitlength,
		uint64_t const rminsplitsize,
		uint64_t const rmaxxdist,
		uint64_t const ractivemax,
		uint64_t const rfracmul,
		uint64_t const rfracdiv,
		uint64_t const ralgndommul,
		uint64_t const ralgndomdiv,
		uint64_t const rchaindommul,
		uint64_t const rchaindomdiv,
		bool const rselfcheck,
		uint64_t const rchainminscore,
		uint64_t const rmaxocc,
		uint64_t const rminalgnlen,
		unsigned int const rmaxwerr,
                int64_t const rmaxback,
		libmaus2::parallel::SynchronousCounter<uint64_t> & rcnt,
		libmaus2::timing::RealTimeClock & rrtc
	) : LW(rLW), PFAI(rPFAI), meta(rmeta), Prank(rPrank), Pcache(rPcache), cocache(Prank,meta), n(Prank.size()), text(rtext),
	    minfreq(rminfreq), minlen(rminlen), limit(rlimit), minsplitlength(rminsplitlength), minsplitsize(rminsplitsize),
	    proc(rmeta,cocache,rPrank,rBSSSA,rtext,rmaxxdist,ractivemax,rfracmul,rfracdiv,rselfcheck,rchainminscore,rmaxocc,ralgndommul,ralgndomdiv,rchaindommul,rchaindomdiv,rmaxwerr,rmaxback),
	    nnp(rmaxwerr,rmaxback),
	    onrefid(0), on(0), off(0), cnt(rcnt), rtc(rrtc),
	    minalgnlen(rminalgnlen)
	{

	}

	void process(
		std::string const & sid,
		std::string const & rpat,
		uint64_t const rid,
		libmaus2::bambam::BamAlignment const * algn
	)
	{
		bool const algnisReverse = algn ? algn->isReverse() : false;
		int64_t const algnrefid = algn ? algn->getRefID() : -1;
		int64_t const algnpos = algn ? algn->getPos() : -1;
		int64_t const algnreflen = algn ? algn->getReferenceLength() : -1;

		std::string pat = libmaus2::fastx::mapString(rpat);

		for ( uint64_t i = 0; i < pat.size(); ++i )
			if ( pat[i] >= 4 )
				pat[i] = libmaus2::random::Random::rand8() % 4;

		uint64_t const Psize = pat.size();
		char const * c = pat.c_str();

		libmaus2::rank::DNARankSMEMComputation::SMEMEnumerator<char const *> senum(
			Prank,
			&Pcache,
			c,
			0,
			Psize,
			minfreq,
			minlen,
			limit,
			minsplitlength,
			minsplitsize
		);

		proc.process(senum,c,Psize);

		// more than one alignment?
		if ( algn && algnrefid >= 0 && (proc.CNIS.aalgno > 1) )
		{
			int64_t beston = -1;
			int64_t firstscore = -1;

			libmaus2::lcs::NNPAlignResult const & res0 = proc.CNIS.Aalgn[0].res;
			libmaus2::fastx::DNAIndexMetaDataBigBandBiDir::Coordinates C0 = meta.mapCoordinatePair(res0.bbpos,res0.bepos);
			libmaus2::math::IntegerInterval<int64_t> I0(C0.left,C0.left+C0.length);
			libmaus2::math::IntegerInterval<int64_t> IR(algnpos,algnpos+algnreflen-1);

			if (
				static_cast<int64_t>(C0.seq) != algnrefid
				||
				static_cast<bool>(C0.rc) != algnisReverse
				||
				I0.intersection(IR).isEmpty()
			)
			{
				libmaus2::parallel::ScopePosixSpinLock slock(libmaus2::aio::StreamLock::cerrlock);
				std::cerr << "[A] " << sid << std::endl;

				if ( static_cast<int64_t>(C0.seq) != algnrefid )
				{
					std::cerr << "[E] wrong refid" << std::endl;
				}
				else if ( C0.rc != algnisReverse )
				{
					std::cerr << "[E] wrong strand" << std::endl;
				}
				else
				{
					std::cerr << "[E] no overlap " << I0 << "\t" << IR << std::endl;
				}

				for ( uint64_t i = 0; i < proc.CNIS.aalgno; ++i )
				{
					libmaus2::lcs::NNPAlignResult const & res = proc.CNIS.Aalgn[i].res;
					libmaus2::fastx::DNAIndexMetaDataBigBandBiDir::Coordinates C = meta.mapCoordinatePair(res.bbpos,res.bepos);
					if ( C.valid )
					{
						bool const seqok = (static_cast<int64_t>(C.seq) == algnrefid);
						bool const rcok = (C.rc == algnisReverse);
						libmaus2::math::IntegerInterval<int64_t> IC(C.left,C.left+C.length);
						bool const Cok = !IC.intersection(IR).isEmpty();
						bool const ok = seqok && rcok && Cok;

						std::cerr << (ok?"*\t":"\t") << res << " " << C.seq << "," << C.left << " actual " << algnrefid << " " << algnpos << " score " << res.getScore() << std::endl;

						if ( ok && beston < 0 )
							beston = res.getScore();
						if ( firstscore < 0 )
							firstscore = res.getScore();
					}
				}
			}

			if ( beston >= 0 && firstscore > 0 )
			{
				libmaus2::parallel::ScopePosixSpinLock slock(libmaus2::aio::StreamLock::cerrlock);
				std::cerr << "best frac " << static_cast<double>(beston) / static_cast<double>(firstscore) << std::endl;
			}
		}

		bool first = true;
		for ( uint64_t i = 0; i < proc.CNIS.aalgno; ++i )
		{
			libmaus2::lcs::NNPAlignResult const & res = proc.CNIS.Aalgn[i].res;
			libmaus2::fastx::DNAIndexMetaDataBigBandBiDir::Coordinates C = meta.mapCoordinatePair(res.bbpos,res.bepos);

			if ( C.valid && (res.aepos-res.abpos >= minalgnlen) )
			{
				std::string const apat = C.rc ? libmaus2::fastx::reverseComplement(pat) : pat;

				uint64_t const seedposref = (C.rc ? (n - proc.CNIS.Aalgn[i].seedposb - proc.CNIS.Aalgn[i].seedlength) : proc.CNIS.Aalgn[i].seedposb) - meta.L[C.seq];
				uint64_t const seedposq = C.rc ? (Psize - proc.CNIS.Aalgn[i].seedposa - proc.CNIS.Aalgn[i].seedlength) : proc.CNIS.Aalgn[i].seedposa;

				char const * seq_ref   = text + meta.L[C.seq];
				char const * seq_query = apat.c_str(); //C.rc ? pat.c_str();

				libmaus2::lcs::NNPAlignResult const nres = nnp.align(
					seq_ref,seq_ref + meta.S[C.seq].l,seedposref,
					seq_query,seq_query+Psize,seedposq,
					trace
				);

				#if 0
				trace.printTraceLines(
					std::cerr,
					seq_ref + nres.abpos,
					seq_query + nres.bbpos,
					80,
					std::string(" "),
					std::string("\n"),
					libmaus2::fastx::remapChar
				);
				#endif

				trace.computeTrace(ATC);

				if ( first )
				{
					uint64_t const ncigar = libmaus2::bambam::CigarStringParser::traceToCigar(
						ATC,Aopblocks,Aop,0,nres.bbpos,Psize-nres.bepos,0
					);

					if ( Psize > Hbuf.size() )
					{
						Hbuf.ensureSize(Psize);
						std::fill(Hbuf.begin(),Hbuf.end(),255);
					}

					libmaus2::bambam::BamAlignmentEncoderBase::encodeAlignmentPreMapped(
						buffer,
						sid.c_str(),
						sid.size(),
						C.seq,
						nres.abpos,
						255 /* mapq */,
						C.rc ? libmaus2::bambam::BamFlagBase::LIBMAUS2_BAMBAM_FREVERSE : 0,
						Aop.begin(),
						ncigar,
						-1 /* next ref id */,
						-1 /* next pos */,
						0 /* tlen */,
						apat.c_str() /* seq */,
						Psize /* seq len */,
						Hbuf.begin() /* qual */,
						0 /* quality offset */
					);

					first = false;
				}
				else
				{
					uint64_t const ncigar = libmaus2::bambam::CigarStringParser::traceToCigar(
						ATC,Aopblocks,Aop,nres.bbpos,0,0,Psize-nres.bepos
					);

					libmaus2::bambam::BamAlignmentEncoderBase::encodeAlignmentPreMapped(
						buffer,
						sid.c_str(),
						sid.size(),
						C.seq,
						nres.abpos,
						255 /* mapq */,
						(C.rc ? libmaus2::bambam::BamFlagBase::LIBMAUS2_BAMBAM_FREVERSE : 0) | libmaus2::bambam::BamFlagBase::LIBMAUS2_BAMBAM_FSECONDARY,
						Aop.begin(),
						ncigar,
						-1 /* next ref id */,
						-1 /* next pos */,
						0 /* tlen */,
						apat.c_str() /* seq */,
						0 /* seq len */,
						apat.c_str() /* qual */
					);
				}

				libmaus2::bambam::BamAlignmentEncoderBase::putAuxNumber(buffer,"AS",'i',res.getScore());
				libmaus2::bambam::BamAlignmentEncoderBase::putAuxNumber(buffer,"NM",'i',res.dif);

				buffer_type::shared_ptr_type Pbuffer = LW.bufferFreeList.get();
				Pbuffer->swap(buffer);
				SQB.push_back(Pbuffer);

				if ( static_cast<int64_t>(C.seq) == algnrefid )
				{
					onrefid++;

					#if 0
					libmaus2::math::IntegerInterval<int64_t> IO(algnpos,algnpos+algnreflen-1);
					libmaus2::math::IntegerInterval<int64_t> IA(C.left,C.left+C.length-1);
					std::cerr << "aligned=" << IA << " refbam=" << IO << std::endl;
					#endif
				}
			}
		}

		// no alignment? output unmapped read
		if ( ! SQB.size() )
		{
			if ( Psize > Hbuf.size() )
			{
				Hbuf.ensureSize(Psize);
				std::fill(Hbuf.begin(),Hbuf.end(),255);
			}

			libmaus2::bambam::BamAlignmentEncoderBase::encodeAlignmentPreMapped(
				buffer,
				sid.c_str(),
				sid.size(),
				-1,
				-1,
				255 /* mapq */,
				libmaus2::bambam::BamFlagBase::LIBMAUS2_BAMBAM_FUNMAP,
				Aop.begin(),
				0 /* ncigar */,
				-1 /* next ref id */,
				-1 /* next pos */,
				0 /* tlen */,
				pat.c_str() /* seq */,
				Psize /* seq len */,
				Hbuf.begin() /* qual */,
				0 /* quality offset */
			);

			buffer_type::shared_ptr_type Pbuffer = LW.bufferFreeList.get();
			Pbuffer->swap(buffer);
			SQB.push_back(Pbuffer);
		}

		LW.put(SQB,rid);

		if ( onrefid )
			on++;
		else
			off++;

		uint64_t const lcnt = ++cnt;
		if ( (lcnt % 256) == 0 )
		{
			double const tim = rtc.getElapsedSeconds();

			libmaus2::parallel::ScopePosixSpinLock slock(libmaus2::aio::StreamLock::cerrlock);

			if ( algn && algnrefid >= 0 )
				std::cerr << lcnt << "\t" << sid << "\t" << proc.CNIS.aalgno << "\t" << onrefid << "\t" << on << "\t" << off
					// << "\t" << algn.getAuxString("er")
					<< "\t"
					<< rtc.formatTime(tim)
					<< "\t"
					<< (lcnt / tim)
					<< std::endl;
			else
				std::cerr
					<< lcnt
					<< "\t" << sid
					<< "\t" << proc.CNIS.aalgno
					<< "\t" << rtc.formatTime(tim)
					<< "\t" << (lcnt / tim)
					<< std::endl;
		}
	}
};


static uint64_t getDefaultNumThreads()
{
	return libmaus2::suffixsort::bwtb3m::BwtMergeSortOptions::getDefaultNumThreads();
}

static std::string getDefaultInputFormat()
{
	return "fasta";
}

static std::string getDefaultOutputFormat()
{
	return "bam";
}

static int64_t getDefaultConstructionMemory()
{
	return 3*libmaus2::util::MemoryStatistics::getPhysicalMemory()/4;
}

static uint64_t getDefaultSuffixArraySamplingRate()
{
	return 32;
}

static uint64_t getDefaultInverseSuffixArraySamplingRate()
{
	return 256;
}

static std::string getDefaultTmpPrefix(std::string const & progname)
{
	return libmaus2::util::ArgInfo::getDefaultTmpFileName(progname);
}

static uint64_t getDefaultMaxOcc()
{
	return 500;
}

static uint64_t getDefaultChainMinScore()
{
	return 20;
}

static uint64_t getDefaultMinFreq()
{
	return 1;
}

static uint64_t getDefaultMinLen()
{
	return 20;
}

static uint64_t getDefaultLimit()
{
	return 32;
}

static uint64_t getDefaultMinSplitLength()
{
	return 28;
}

static uint64_t getDefaultMinSplitSize()
{
	return 10;
}

static uint64_t getDefaultMaxXDist()
{
	return 1000;
}

static uint64_t getDefaultActiveMax()
{
	return 1;
}

static uint64_t getDefaultFracMul()
{
	return 95;
}

static uint64_t getDefaultFracDiv()
{
	return 100;
}

static uint64_t getDefaultChainDomMul()
{
	return 95;
}

static uint64_t getDefaultChainDomDiv()
{
	return 100;
}

static uint64_t getDefaultAlignDomMul()
{
	return 95;
}

static uint64_t getDefaultAlignDomDiv()
{
	return 100;
}

static uint64_t getDefaultMinAlgnLen()
{
	return 50;
}

static uint64_t getDefaultCacheK()
{
	return 12;
}

#if 0
static std::string formatNumber(int64_t const v)
{
	std::ostringstream ostr;
	ostr << v;
	return ostr.str();
}
#endif

template<typename default_type>
static std::string formatRHS(std::string const & description, default_type def)
{
	std::ostringstream ostr;
	ostr << description << " (default " << def << ")";
	return ostr.str();
}

/*
 parameters:

 -i : default fasta, inputformat
 -o : default bam, outputformat
 -t : default number of logical cores, threads
 --constructmem: default 16GB, BWT+SA construction memory
 --sasamplingrate: default 32, sampling rate for suffix array
 --isasamplingrate: default 256, sampling rate for inverse suffix array
 -T : tmp prefix
 --maxocc: default 50, frequency threshold for seeds (seeds more frequent will be sub sampled)
 --chainminscore: default 20, minimum chain score (chains with smaller score are discarded)
 --minfreq: default 1, minimum frequency for seeds
 --minlen: default 20, minimum seed length
 --limit: default 32, maximum seed extension left and right
 --minsplitlength: default 28, minimum seed length to consider splitting
 --minsplitsize: default 10, maximum frequency to consider splitting
 --maxxdist: default 1000, maximum distance to consider seeds chainable
 --activemax: default 1, minimum number of chains covering another chain to consider it dominated
 --fracmul: default 95, fragment overlap threshold counter
 --fracdiv: default 100, fragment overlap threshold denominator
 --chaindommul: default 95, chain domination threshold counter
 --chaindomdiv: default 100, chain domination threshold denominator
 --algndommul: default 95, alignment domination threshold counter
 --algndomdiv: default 100, alignment domination threshold denominator
 --minalgnlen: default 50, minimum length of alignment reported (length in reference bases)
 -K: kmer cache K size
 */

static std::string helpMessage(libmaus2::util::ArgParser const & arg)
{
	std::vector < std::pair < std::string, std::string > > optionMap;
	optionMap . push_back ( std::pair < std::string, std::string >("i", formatRHS("inputformat",getDefaultInputFormat())));
	optionMap . push_back ( std::pair < std::string, std::string >("o", formatRHS("inputformat",getDefaultOutputFormat())));
	optionMap . push_back ( std::pair < std::string, std::string >("t", formatRHS("number of threads",getDefaultNumThreads())));
	optionMap . push_back ( std::pair < std::string, std::string >("constructmem", formatRHS("memory guide for BWT+SA construction",getDefaultConstructionMemory())));
	optionMap . push_back ( std::pair < std::string, std::string >("sasamplingrate", formatRHS("suffix array sampling rate",getDefaultSuffixArraySamplingRate())));
	optionMap . push_back ( std::pair < std::string, std::string >("isasamplingrate", formatRHS("inverse suffix array sampling rate",getDefaultInverseSuffixArraySamplingRate())));
	optionMap . push_back ( std::pair < std::string, std::string >("T", formatRHS("temporary file prefix",getDefaultTmpPrefix(arg.progname))));
	optionMap . push_back ( std::pair < std::string, std::string >("minalgnlen", formatRHS("minimum length of alignments reported",getDefaultMinAlgnLen())));
	optionMap . push_back ( std::pair < std::string, std::string >("maxocc", formatRHS("frequency threshold for seeds",getDefaultMaxOcc())));
	optionMap . push_back ( std::pair < std::string, std::string >("chainminscore", formatRHS("minimum chain score",getDefaultChainMinScore())));
	optionMap . push_back ( std::pair < std::string, std::string >("minfreq", formatRHS("minimum frequency for seeds",getDefaultMinFreq())));
	optionMap . push_back ( std::pair < std::string, std::string >("minlen", formatRHS("minimum seed length",getDefaultMinLen())));
	optionMap . push_back ( std::pair < std::string, std::string >("limit", formatRHS("maximum seed extension left and right",getDefaultLimit())));
	optionMap . push_back ( std::pair < std::string, std::string >("minsplitlength", formatRHS("minimum seed length to consider splitting",getDefaultMinSplitLength())));
	optionMap . push_back ( std::pair < std::string, std::string >("minsplitsize", formatRHS("maximum frequency to consider splitting",getDefaultMinSplitSize())));
	optionMap . push_back ( std::pair < std::string, std::string >("maxxdist", formatRHS("maximum distance to consider seeds chainable",getDefaultMaxXDist())));
	optionMap . push_back ( std::pair < std::string, std::string >("activemax", formatRHS("minimum number of chains covering another chain to consider it dominated",getDefaultActiveMax())));
	optionMap . push_back ( std::pair < std::string, std::string >("fracmul", formatRHS("fragment overlap threshold counter",getDefaultFracMul())));
	optionMap . push_back ( std::pair < std::string, std::string >("fracdiv", formatRHS("fragment overlap threshold denominator",getDefaultFracDiv())));
	optionMap . push_back ( std::pair < std::string, std::string >("chaindommul", formatRHS("chain domination threshold counter",getDefaultChainDomMul())));
	optionMap . push_back ( std::pair < std::string, std::string >("chaindomdiv", formatRHS("chain domination threshold denominator",getDefaultChainDomDiv())));
	optionMap . push_back ( std::pair < std::string, std::string >("algndommul", formatRHS("alignment domination threshold counter",getDefaultAlignDomMul())));
	optionMap . push_back ( std::pair < std::string, std::string >("algndomdiv", formatRHS("alignment domination threshold denominator",getDefaultAlignDomDiv())));
	optionMap . push_back ( std::pair < std::string, std::string >("K", formatRHS("kmer cache K size",getDefaultCacheK())));

	uint64_t maxlhs = 0;
	for ( std::vector < std::pair < std::string, std::string > >::const_iterator ita = optionMap.begin(); ita != optionMap.end(); ++ita )
	{
		assert ( ita->first.size() );

		if ( ita->first.size() == 1 )
			maxlhs = std::max(maxlhs,static_cast<uint64_t>(ita->first.size()+1));
		else
			maxlhs = std::max(maxlhs,static_cast<uint64_t>(ita->first.size()+2));
	}

	std::ostringstream messtr;
	for ( std::vector < std::pair < std::string, std::string > >::const_iterator ita = optionMap.begin(); ita != optionMap.end(); ++ita )
	{
		std::string const key = ita->first;

		messtr << "\t";
		messtr << std::setw(maxlhs) << std::setfill(' ');
		if ( key.size() == 1 )
			messtr << (std::string("-")+key);
		else
			messtr << (std::string("--")+key);

		messtr << std::setw(0);

		messtr << ": ";

		messtr << ita->second;
		messtr << "\n";
	}

	return messtr.str();
}

int yalora(libmaus2::util::ArgParser const & arg, std::string const & fn)
{
	// number of threas
	uint64_t const numthreads = arg.uniqueArgPresent("t") ? arg.getUnsignedNumericArg<uint64_t>("t") : getDefaultNumThreads();
	// input format
	std::string const inputformat = arg.uniqueArgPresent("i") ? arg["i"] : getDefaultInputFormat();
	// output format
	std::string const outputformat = arg.uniqueArgPresent("o") ? arg["o"] : getDefaultOutputFormat();

	// memory used for bwt construction
	uint64_t const constructmem = arg.uniqueArgPresent("constructmem") ? arg.getUnsignedNumericArg<uint64_t>("constructmem") : getDefaultConstructionMemory();
	// SA sampling rate
	uint64_t const sasamplingrate = arg.uniqueArgPresent("sasamplingrate") ? arg.getUnsignedNumericArg<uint64_t>("sasamplingrate") : getDefaultSuffixArraySamplingRate();
	// ISA sampling rate
	uint64_t const isasamplingrate = arg.uniqueArgPresent("isasamplingrate") ? arg.getUnsignedNumericArg<uint64_t>("isasamplingrate") : getDefaultInverseSuffixArraySamplingRate();

	// tmp prefix format
	std::string const tmpprefix = arg.uniqueArgPresent("T") ? arg["T"] : getDefaultTmpPrefix(arg.progname);

	bool const selfcheck = false;
	//uint64_t const chainminscore = 20;
	uint64_t const maxocc = arg.uniqueArgPresent("maxocc") ? arg.getUnsignedNumericArg<uint64_t>("maxocc") : getDefaultMaxOcc();
	uint64_t const chainminscore = arg.uniqueArgPresent("chainminscore") ? arg.getUnsignedNumericArg<uint64_t>("chainminscore") : getDefaultChainMinScore();

	// minimum frequency for seeds
	uint64_t const minfreq = arg.uniqueArgPresent("minfreq") ? arg.getUnsignedNumericArg<uint64_t>("minfreq") : getDefaultMinFreq();

	// minimum seed length
	uint64_t const minlen = arg.uniqueArgPresent("minlen") ? arg.getUnsignedNumericArg<uint64_t>("minlen") : getDefaultMinLen();

	// maximum extension to left and right for SMEMs
	uint64_t const limit = arg.uniqueArgPresent("limit") ? arg.getUnsignedNumericArg<uint64_t>("limit") : getDefaultLimit();

	// seed splitting
	uint64_t const minsplitlength = arg.uniqueArgPresent("minsplitlength") ? arg.getUnsignedNumericArg<uint64_t>("minsplitlength") : getDefaultMinSplitLength();
	uint64_t const minsplitsize = arg.uniqueArgPresent("minsplitsize") ? arg.getUnsignedNumericArg<uint64_t>("minsplitsize") : getDefaultMinSplitSize();

	// maximum distance in X/Y to connect up a seed to a previously existing chain
	uint64_t const maxxdist = arg.uniqueArgPresent("maxxdist") ? arg.getUnsignedNumericArg<uint64_t>("maxxdist") : getDefaultMaxXDist();
	//
	uint64_t const activemax = arg.uniqueArgPresent("activemax") ? arg.getUnsignedNumericArg<uint64_t>("activemax") : getDefaultActiveMax();

	// overlap fraction
	uint64_t const fracmul = arg.uniqueArgPresent("fracmul") ? arg.getUnsignedNumericArg<uint64_t>("fracmul") : getDefaultFracMul();
	uint64_t const fracdiv = arg.uniqueArgPresent("fracdiv") ? arg.getUnsignedNumericArg<uint64_t>("fracdiv") : getDefaultFracDiv();

	// chain dominance fraction
	uint64_t const chaindommul = arg.uniqueArgPresent("chaindommul") ? arg.getUnsignedNumericArg<uint64_t>("chaindommul") : getDefaultChainDomMul();
	uint64_t const chaindomdiv = arg.uniqueArgPresent("chaindomdiv") ? arg.getUnsignedNumericArg<uint64_t>("chaindomdiv") : getDefaultChainDomDiv();

	// alignment dominance fraction
	uint64_t const algndommul = arg.uniqueArgPresent("algndommul") ? arg.getUnsignedNumericArg<uint64_t>("algndommul") : getDefaultAlignDomMul();
	uint64_t const algndomdiv = arg.uniqueArgPresent("algndomdiv") ? arg.getUnsignedNumericArg<uint64_t>("algndomdiv") : getDefaultAlignDomDiv();

	// minimum length of stored alignment on reference
	uint64_t const minalgnlen = arg.uniqueArgPresent("minalgnlen") ? arg.getUnsignedNumericArg<uint64_t>("minalgnlen") : getDefaultMinAlgnLen();

	// kmer cache K
	uint64_t const cachek = arg.uniqueArgPresent("K") ? arg.getUnsignedNumericArg<uint64_t>("K") : getDefaultCacheK();

	// maximum number of errors in 64 symbol window
	unsigned int const maxwerr =
		arg.uniqueArgPresent("maxwerr") ? arg.getUnsignedNumericArg<uint64_t>("maxwerr") :
			libmaus2::lcs::NNP::getDefaultMaxWindowError();
	// distance threshold from best trace
	int64_t const maxback =
		arg.uniqueArgPresent("maxback") ? arg.getUnsignedNumericArg<uint64_t>("maxback") :
			libmaus2::lcs::NNP::getDefaultMaxBack();

	// set up arguments for constructing output writer
	libmaus2::util::ArgInfo warginfo(arg.progname);
	warginfo.insertKey("outputformat",outputformat);
	warginfo.insertKey("level","0");
	warginfo.insertKey("reference",fn);

	std::string const compactfn = fn + ".compact";
	std::string const compactmetafn = compactfn + ".meta";

	// produce compact file if it does not exist
	if ( ! libmaus2::util::GetFileSize::fileExists(compactfn) || libmaus2::util::GetFileSize::isOlder(compactfn,fn) )
	{
		libmaus2::fastx::FastAToCompact4BigBandBiDir::fastaToCompact4BigBandBiDir(
			std::vector<std::string>(1,fn),
			&(std::cerr),
			false /* single strand */,
			compactfn
		);
	}

	// bwt name
	std::string const bwtfn = fn + ".bwt";
	// bwt stats name
	std::string const bwtmetafn = bwtfn + ".meta";

	// construct BWT if it does not exist, otherwise load info
	libmaus2::suffixsort::bwtb3m::BwtMergeSortResult res;
	if ( ! libmaus2::util::GetFileSize::fileExists(bwtmetafn) || libmaus2::util::GetFileSize::isOlder(bwtmetafn,compactfn) )
	{
		libmaus2::suffixsort::bwtb3m::BwtMergeSortOptions options(
			compactfn,
			constructmem, // mem
			// libmaus2::suffixsort::bwtb3m::BwtMergeSortOptions::getDefaultMem(),
			numthreads,
			"compactstream",
			false /* bwtonly */,
			std::string(tmpprefix+"_bwtb3m"),
			std::string(), // sparse
			bwtfn,
			isasamplingrate /* isa */,
			sasamplingrate /* sa */
		);

		res = libmaus2::suffixsort::bwtb3m::BwtMergeSort::computeBwt(options,&std::cerr);
		res.serialise(bwtmetafn);
	}
	else
	{
		res.deserialise(bwtmetafn);
	}

	// generate FastA index if it does not exist
	std::string const fainame = fn+".fai";
	libmaus2::fastx::FastAIndexGenerator::generate(fn,fainame,true /* verbose */);

	// generate sequence md5 file if it does not exist
	std::string const m5info = fn+".m5info";
	if ( ! libmaus2::util::GetFileSize::fileExists(m5info) || libmaus2::util::GetFileSize::isOlder(m5info,fn) )
	{
		libmaus2::aio::InputStreamInstance ISI(fn);
		libmaus2::fastx::FastAStreamSet FASS(ISI);
		std::map<std::string,std::string> m5map = FASS.computeMD5(false,false);

		libmaus2::aio::OutputStreamInstance OSI(m5info);
		for ( std::map<std::string,std::string>::const_iterator ita = m5map.begin(); ita != m5map.end(); ++ita )
		{
			libmaus2::util::StringSerialisation::serialiseString(OSI,ita->first);
			libmaus2::util::StringSerialisation::serialiseString(OSI,ita->second);
		}
	}

	// load sequence MD5 file
	std::map<std::string,std::string> m5map;
	{
		libmaus2::aio::InputStreamInstance ISI(m5info);
		while ( ISI && ISI.peek() != std::istream::traits_type::eof() )
		{
			std::string const key = libmaus2::util::StringSerialisation::deserialiseString(ISI);
			std::string const value = libmaus2::util::StringSerialisation::deserialiseString(ISI);
			m5map[key] = value;
		}
	}

	// load FastA index
	libmaus2::fastx::FastAIndex::unique_ptr_type PFAI(libmaus2::fastx::FastAIndex::load(fainame));

	// load index meta data
	libmaus2::fastx::DNAIndexMetaDataBigBandBiDir::unique_ptr_type Pmeta(libmaus2::fastx::DNAIndexMetaDataBigBandBiDir::load(compactmetafn));
	std::ostringstream headerostr;
	headerostr << "@HD\tVN:1.5\tSO:unknown\n";
	for ( uint64_t i = 0; i < Pmeta->S.size()/2; ++i )
	{
		headerostr << "@SQ\tLN:" << Pmeta->S[i].l << "\tSN:" << (*PFAI)[i].name;
		if ( m5map.find((*PFAI)[i].name) != m5map.end() )
			headerostr << "\tM5:" << m5map.find((*PFAI)[i].name)->second;
		headerostr << "\n";
	}
	headerostr << "@PG\tID:yalora\tPN:yalora\tCL:" << arg.commandline << "\tVN:" << PACKAGE_VERSION << "\n";
	// deallocate m5 map
	m5map.clear();
	// construct BAM header
	libmaus2::bambam::BamHeader header(headerostr.str());

	// open output
	libmaus2::bambam::BamBlockWriterBase::unique_ptr_type Pwriter(libmaus2::bambam::BamBlockWriterBaseFactory::construct(header,warginfo));
	libmaus2::bambam::BamBlockWriterBase & writer = *Pwriter;

	// reorder buffer for output
	LockedWriter LW(writer,tmpprefix+"_output_tmp");

	// load index (BWT+sampled SA)
	libmaus2::rank::DNARank::unique_ptr_type Prank(res.loadDNARank(numthreads));
	// length of string
	uint64_t const n = Prank->size();
	#if 0
	typedef libmaus2::suffixsort::bwtb3m::BwtMergeSortResult::BareSimpleSampledSuffixArray ssa_type;
	ssa_type BSSSA(res.loadBareSimpleSuffixArray());
	#else
	typedef libmaus2::suffixsort::bwtb3m::BwtMergeSortResult::CompactBareSimpleSampledSuffixArray ssa_type;
	ssa_type BSSSA(res.loadCompactBareSimpleSuffixArray(n));
	#endif
	// load coordinate cache
	libmaus2::fastx::CoordinateCacheBiDir cocache(*Prank,*Pmeta,16 /* blockshfit */);

	// load text
	libmaus2::autoarray::AutoArray<char> A(n,false);
	{
		libmaus2::bitio::CompactDecoderWrapper CDW(compactfn);
		CDW.read(A.begin(),n);
		assert ( CDW.gcount() == static_cast<int64_t>(n) );
	}

	libmaus2::rank::DNARankKmerCache::unique_ptr_type Pcache(new libmaus2::rank::DNARankKmerCache(*Prank,cachek,numthreads));

	libmaus2::parallel::SynchronousCounter<uint64_t> cnt;
	libmaus2::autoarray::AutoArray < AlignContext<ssa_type>::unique_ptr_type > Acontext(numthreads);
	libmaus2::timing::RealTimeClock rtc;
	rtc.start();

	for ( uint64_t i = 0; i < numthreads; ++i )
	{
		AlignContext<ssa_type>::unique_ptr_type tcontext(new AlignContext<ssa_type>(LW,*PFAI,*Pmeta,*Prank,*Pcache,BSSSA,A.begin(),minfreq,minlen,limit,minsplitlength,minsplitsize,maxxdist,activemax,fracmul,fracdiv,algndommul,algndomdiv,chaindommul,chaindomdiv,selfcheck,chainminscore,maxocc,minalgnlen,maxwerr,maxback,cnt,rtc));
		Acontext[i] = UNIQUE_PTR_MOVE(tcontext);
	}


	LockedGetInterface::unique_ptr_type PLG;

	if ( inputformat == "bam" )
	{
		LockedGetInterface::unique_ptr_type TLG(new BamLockedGet(std::cin));
		PLG = UNIQUE_PTR_MOVE(TLG);
	}
	else
	{
		LockedGetInterface::unique_ptr_type TLG(new FastALockedGet(std::cin));
		PLG = UNIQUE_PTR_MOVE(TLG);
	}

	int volatile ret = EXIT_SUCCESS;

	LockedGetInterface & LG = *PLG;

	#if defined(_OPENMP)
	#pragma omp parallel num_threads(numthreads)
	#endif
	{
		uint64_t const tid =
			#if defined(_OPENMP)
			omp_get_thread_num()
			#else
			0
			#endif
			;

		libmaus2::bambam::BamAlignment algn;
		libmaus2::bambam::BamAlignment * palgn = 0;
		uint64_t rid;
		std::string sid;
		std::string rpat;

		try
		{
			while ( LG.getNext(sid,rpat,rid,algn,palgn) )
			{
				// Acontext[tid]->process(sid,rpat,rid,&algn);
				Acontext[tid]->process(sid,rpat,rid,palgn);
			}
		}
		catch(std::exception const & ex)
		{
			libmaus2::parallel::ScopePosixSpinLock slock(libmaus2::aio::StreamLock::cerrlock);
			std::cerr << ex.what() << std::endl;

			ret = EXIT_FAILURE;
		}
	}

	LW.flush();

	return ret;
}

int main(int argc, char * argv[])
{
	try
	{
		libmaus2::util::ArgParser arg(argc,argv);

		if ( arg.uniqueArgPresent("v") || arg.uniqueArgPresent("version") )
		{
			std::cerr << "This is " << PACKAGE_NAME << " version " << PACKAGE_VERSION << "." << std::endl;
			std::cerr << PACKAGE_NAME << " is distributed under version 3 of the GPL." << std::endl;
			return EXIT_SUCCESS;
		}
		else if ( arg.uniqueArgPresent("h") || arg.uniqueArgPresent("help") || arg.size() < 1 )
		{
			std::cerr << "This is " << PACKAGE_NAME << " version " << PACKAGE_VERSION << "." << std::endl;
			std::cerr << PACKAGE_NAME << " is distributed under version 3 of the GPL." << std::endl;
			std::cerr << "\n";
			std::cerr << "usage: " << arg.progname << " [options] ref.fasta < input.fasta\n";
			std::cerr << "\n";
			std::cerr << "The following options can be used (no space between option name and parameter allowed):\n\n";
			std::cerr << helpMessage(arg);
			return EXIT_SUCCESS;
		}
		else
		{
			return yalora(arg,arg[0]);
		}
	}
	catch(std::exception const & ex)
	{
		std::cerr << ex.what() << std::endl;
		return EXIT_FAILURE;
	}
}
