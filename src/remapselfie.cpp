/*
    yalora
    Copyright (C) 2017 German Tischler

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
#include <config.h>

#include <Coord.hpp>
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
#include <libmaus2/lcs/NNPCor.hpp>
#include <libmaus2/lcs/NNP.hpp>
#include <libmaus2/lcs/NP.hpp>
#include <libmaus2/lcs/NPLLinMem.hpp>
#include <libmaus2/lcs/AlignmentPrint.hpp>
#include <libmaus2/dazzler/align/Overlap.hpp>
#include <libmaus2/dazzler/align/AlignmentWriter.hpp>
#include <libmaus2/dazzler/align/AlignmentWriterArray.hpp>
#include <libmaus2/dazzler/align/AlignmentFile.hpp>


static bool getLine(std::istream & in, std::string & line)
{
	std::getline(in,line);

	return in && line.size();
}

struct SortEntry
{
	// absolute position in index
	uint64_t zzz;
	// sequence
	uint64_t seq;
	// position in seq
	uint64_t pos;
	// length of hit
	uint64_t l;
	// is entry active?
	bool active;

	SortEntry()
	{

	}
	SortEntry(uint64_t const rzzz, uint64_t const rseq, uint64_t const rpos, uint64_t const rl)
	: zzz(rzzz), seq(rseq), pos(rpos), l(rl), active(true) {}

	int64_t getDiagonal() const
	{
		int64_t izzz = zzz;
		int64_t ipos = pos;
		return izzz - ipos;
	}

	bool overlap(SortEntry const & S) const
	{
		if ( getDiagonal() != S.getDiagonal() )
			return false;

		libmaus2::math::IntegerInterval<int64_t> IA(  zzz,  zzz+  l-1);
		libmaus2::math::IntegerInterval<int64_t> IB(S.zzz,S.zzz+S.l-1);

		return !IA.intersection(IB).isEmpty();
	}

	int64_t getAntiDiagonalLow() const
	{
		return (zzz + pos);
	}

	int64_t getAntiDiagonalHigh() const
	{
		return (zzz + pos + 2*l);
	}

	libmaus2::math::IntegerInterval<int64_t> getAntiDiagonalInterval() const
	{
		return libmaus2::math::IntegerInterval<int64_t>(
			getAntiDiagonalLow(),
			getAntiDiagonalHigh()-1
		);
	}

	bool operator<(SortEntry const & S) const
	{
		if ( seq != S.seq )
			return seq < S.seq;

		int64_t const d0 = getDiagonal();
		int64_t const d1 = S.getDiagonal();

		if ( d0 != d1 )
			return d0 < d1;

		return zzz < S.zzz;
	}

	int64_t zzzendpoint() const
	{
		return zzz + l;
	}

	int64_t posendpoint() const
	{
		return pos + l;
	}
};

std::ostream & operator<<(std::ostream & out, SortEntry const & S)
{
	out << "SortEntry(" << S.zzz << "," << S.seq << "," << S.pos << "," << S.l << ")";
	return out;
}

struct SortEntryZZZComparator
{
	bool operator()(SortEntry const & A, SortEntry const & B) const
	{
		return A.zzz < B.zzz;
	}
};

struct SortEntryPosComparator
{
	bool operator()(SortEntry const & A, SortEntry const & B) const
	{
		return A.pos < B.pos;
	}
};

struct SortEntryDiagComparator
{
	bool operator()(SortEntry const & A, SortEntry const & B) const
	{
		return
			A.getDiagonal()
			<
			B.getDiagonal()
			;
	}
};

void gapFilterZZZ(SortEntry * Asort, uint64_t const ilow, uint64_t const ihigh, int64_t const maxgap, std::vector < std::pair<uint64_t,uint64_t> > & R)
{
	std::sort(
		Asort + ilow,
		Asort + ihigh,
		SortEntryZZZComparator()
	);

	uint64_t gzzzlow = ilow;
	while ( gzzzlow < ihigh )
	{
		uint64_t gzzzhigh = gzzzlow + 1;

		while ( gzzzhigh < ihigh && static_cast<int64_t>(Asort[gzzzhigh].zzz) - Asort[gzzzhigh-1].zzzendpoint() <= maxgap )
			++gzzzhigh;

		R.push_back(std::pair<uint64_t,uint64_t>(gzzzlow,gzzzhigh));

		gzzzlow = gzzzhigh;
	}
}

void gapFilterPos(SortEntry * Asort, uint64_t const ilow, uint64_t const ihigh, int64_t const maxgap, std::vector < std::pair<uint64_t,uint64_t> > & R)
{
	std::sort(
		Asort + ilow,
		Asort + ihigh,
		SortEntryPosComparator()
	);

	uint64_t gzzzlow = ilow;
	while ( gzzzlow < ihigh )
	{
		uint64_t gzzzhigh = gzzzlow + 1;

		while ( gzzzhigh < ihigh && static_cast<int64_t>(Asort[gzzzhigh].pos) - Asort[gzzzhigh-1].posendpoint() <= maxgap )
			++gzzzhigh;

		R.push_back(std::pair<uint64_t,uint64_t>(gzzzlow,gzzzhigh));

		gzzzlow = gzzzhigh;
	}
}

std::vector < std::pair<uint64_t,uint64_t> > gapFilter(SortEntry * Asort, uint64_t const ilow, uint64_t const ihigh, int64_t const maxgap)
{
	std::vector < std::pair<uint64_t,uint64_t> > R;
	R.push_back(std::pair<uint64_t,uint64_t>(ilow,ihigh));

	bool changed = true;

	while ( changed )
	{
		changed = false;

		std::vector < std::pair<uint64_t,uint64_t> > RZZZ;
		for ( uint64_t i = 0; i < R.size(); ++i )
			gapFilterZZZ(Asort,R[i].first,R[i].second,maxgap,RZZZ);

		std::vector < std::pair<uint64_t,uint64_t> > Rpos;
		for ( uint64_t i = 0; i < R.size(); ++i )
			gapFilterPos(Asort,RZZZ[i].first,RZZZ[i].second,maxgap,Rpos);

		changed = Rpos != R;

		R = Rpos;
	}

	return R;
}

static std::string getDefaultTmpPrefix(std::string const & progname)
{
	return libmaus2::util::ArgInfo::getDefaultTmpFileName(progname);
}

static void extendRight(
	char const * a,
	int64_t const an,
	char const * b,
	int64_t const bn,
	libmaus2::lcs::NNPAlignResult & algn,
	libmaus2::lcs::AlignmentTraceContainer & ATC,
	int64_t const ext = 256,
	double const e = 0.4
)
{
	libmaus2::lcs::NPLLinMem npl;
	bool running = true;

	while ( running )
	{
		running = false;

		int64_t const ava = an - algn.aepos;
		int64_t const avb = bn - algn.bepos;

		int64_t const usea = std::min(ava,ext);
		int64_t const useb = std::min(avb,ext);

		npl.np(
			a+algn.aepos,a+algn.aepos+usea,
			b+algn.bepos,b+algn.bepos+useb
		);

		npl.lastGoodWindowBack(64,e);

		std::pair<int64_t,int64_t> const SL = npl.getStringLengthUsed();

		if ( SL.first + SL.second )
		{
			if ( algn.aepos+SL.first != algn.bepos+SL.second )
			{
				#if 0
				if ( SL.first >= ext/2 && SL.second >= ext/2 )
				{
					std::pair<int64_t,int64_t> const adv = npl.advanceMaxA(ext/2);
					std::pair<int64_t,int64_t> const SLex = npl.getStringLengthUsed(npl.ta,npl.ta+adv.second);

					npl.te = npl.ta + adv.second;

					std::cerr << "extend right(1) " << npl.getAlignmentStatistics() << std::endl;

					ATC.push(npl);
					algn.aepos += SLex.first;
					algn.bepos += SLex.second;
				}
				else
				#endif
				{
					// std::cerr << "extend right(2) " << npl.getAlignmentStatistics() << std::endl;

					ATC.push(npl);
					algn.aepos += SL.first;
					algn.bepos += SL.second;
				}

				running = true;
			}
		}
	}
}

static void extendLeft(
	char const * a,
	char const * b,
	libmaus2::lcs::NNPAlignResult & algn,
	libmaus2::lcs::AlignmentTraceContainer & ATC,
	int64_t const ext = 256,
	double const e = 0.4
)
{
	libmaus2::lcs::NPLLinMem npl;
	bool running = true;

	while ( running )
	{
		running = false;

		int64_t const ava = algn.abpos;
		int64_t const avb = algn.bbpos;

		int64_t const usea = std::min(ava,ext);
		int64_t const useb = std::min(avb,ext);

		std::reverse_iterator<char const *> ra(a + algn.abpos);
		std::reverse_iterator<char const *> rae(a + algn.abpos - usea);
		std::reverse_iterator<char const *> rb(b + algn.bbpos);
		std::reverse_iterator<char const *> rbe(b + algn.bbpos - useb);

		npl.np(ra,rae,rb,rbe);

		npl.lastGoodWindowBack(64,e);

		std::pair<int64_t,int64_t> const SL = npl.getStringLengthUsed();

		if ( SL.first + SL.second )
		{
			if ( algn.abpos-SL.first != algn.bbpos-SL.second )
			{
				#if 0
				if ( SL.first >= ext/2 && SL.second >= ext/2 )
				{
					std::pair<int64_t,int64_t> const adv = npl.advanceMaxA(ext/2);
					std::pair<int64_t,int64_t> const SLex = npl.getStringLengthUsed(
						npl.ta,
						npl.ta+adv.second
					);

					npl.te = npl.ta + adv.second;

					std::reverse(npl.ta,npl.te);

					std::cerr << "extend left(1) " << npl.getAlignmentStatistics() << std::endl;

					ATC.prepend(npl.ta,npl.te);

					algn.abpos -= SLex.first;
					algn.bbpos -= SLex.second;
				}
				else
				#endif
				{
					std::reverse(npl.ta,npl.te);
					// std::cerr << "extend left(2) " << npl.getAlignmentStatistics() << std::endl;

					ATC.prepend(npl.ta,npl.te);
					algn.abpos -= SL.first;
					algn.bbpos -= SL.second;
				}

				running = true;
			}
		}
	}
}

static void extend(
	char const * a,
	int64_t const an,
	char const * b,
	int64_t const bn,
	libmaus2::lcs::NNPAlignResult & algn,
	libmaus2::lcs::AlignmentTraceContainer & ATC,
	int64_t const ext = 256,
	double const e = 0.4
)
{
	extendLeft(a,b,algn,ATC,ext,e);
	extendRight(a,an,b,bn,algn,ATC,ext,e);
}

int remapselfie(libmaus2::util::ArgParser const & arg)
{
	std::string const outfn = arg[0];
	std::string const fn = arg[1];
	std::string const selfie = arg[2];

	std::string const compactfn = fn + ".compact";
	std::string const compactmetafn = compactfn + ".meta";

	if (
		! libmaus2::util::GetFileSize::fileExists(compactfn)
		||
		libmaus2::util::GetFileSize::isOlder(compactfn,fn)
	)
	{
		libmaus2::fastx::FastAToCompact4BigBandBiDir::fastaToCompact4BigBandBiDir(
			std::vector<std::string>(1,fn),
			&(std::cerr),
			false /* single strand */,
			compactfn
		);
	}

	std::string const tmpprefix = arg.uniqueArgPresent("T") ? arg["T"] : getDefaultTmpPrefix(arg.progname);

	uint64_t const numthreads =
		arg.uniqueArgPresent("t") ? arg.getUnsignedNumericArg<uint64_t>("t") : libmaus2::suffixsort::bwtb3m::BwtMergeSortOptions::getDefaultNumThreads();

	std::string const bwtfn = fn + ".bwt";
	std::string const bwtmetafn = bwtfn + ".meta";

	libmaus2::suffixsort::bwtb3m::BwtMergeSortResult res;
	if (
		! libmaus2::util::GetFileSize::fileExists(bwtmetafn)
		||
		libmaus2::util::GetFileSize::isOlder(bwtmetafn,compactfn)
	)
	{
		libmaus2::suffixsort::bwtb3m::BwtMergeSortOptions options(
			compactfn,
			16*1024ull*1024ull*1024ull, // mem
			// libmaus2::suffixsort::bwtb3m::BwtMergeSortOptions::getDefaultMem(),
			numthreads,
			"compactstream",
			false /* bwtonly */,
			std::string("mem:tmp_"),
			std::string(), // sparse
			bwtfn,
			16 /* isa */,
			16 /* sa */
		);

		res = libmaus2::suffixsort::bwtb3m::BwtMergeSort::computeBwt(options,&std::cerr);
		res.serialise(bwtmetafn);
	}
	else
	{
		res.deserialise(bwtmetafn);
	}

	//libmaus2::fastx::FastAIndex::unique_ptr_type PFAI(libmaus2::fastx::FastAIndex::load(fn+".fai"));
	libmaus2::fastx::DNAIndexMetaDataBigBandBiDir::unique_ptr_type Pmeta(libmaus2::fastx::DNAIndexMetaDataBigBandBiDir::load(compactmetafn));

	libmaus2::rank::DNARank::unique_ptr_type Prank(res.loadDNARank(numthreads));
	libmaus2::suffixsort::bwtb3m::BwtMergeSortResult::BareSimpleSampledSuffixArray BSSSA(res.loadBareSimpleSuffixArray());
	libmaus2::rank::DNARankGetPosition< libmaus2::suffixsort::bwtb3m::BwtMergeSortResult::BareSimpleSampledSuffixArray > GP(*Prank,BSSSA);
	libmaus2::fastx::CoordinateCacheBiDir cocache(*Prank,*Pmeta,16 /* blockshfit */);
	uint64_t const cachek = 14;
	uint64_t const k = cachek + 2;
	libmaus2::rank::DNARankKmerCache kcache(*Prank,cachek,numthreads);

	uint64_t const numseq = Pmeta->S.size();

	uint64_t const n = Prank->size();
	libmaus2::autoarray::AutoArray<char> A(n,false);
	libmaus2::bitio::CompactDecoderWrapper CDW(compactfn);
	CDW.read(A.begin(),n);
	assert ( CDW.gcount() == static_cast<int64_t>(n) );

	uint64_t const minreport = 100;
	int64_t const tspace = libmaus2::dazzler::align::AlignmentFile::getMinimumNonSmallTspace();
	int64_t const maxgap = 1000;

	libmaus2::aio::InputStreamInstance ISI(selfie);

	std::string line;

	std::vector < Coord > V;

	while ( getLine(ISI,line) )
	{
		Coord coord;

		if ( coord.parse(line) )
		{

			if ( coord.rc )
				continue;

			//std::cerr << coord << std::endl;

			if ( ! V.size() || !V.back().overlap(coord) )
				V.push_back(coord);
			else
				V.back() = V.back().merge(coord);
		}
		else
		{
			std::cerr << "Cannot parse " << line << std::endl;
		}
	}

	std::cerr << "numseq=" << numseq << std::endl;

	for ( uint64_t i = 0; i < numseq/2; ++i )
	{
		char const * a = A.begin() + Pmeta->L[i];
		char const * b = A.begin() + Pmeta->L[numseq-i];
		uint64_t const l = Pmeta->L[i+1]-Pmeta->L[i];

		for ( uint64_t j = 0; j < l; ++j )
		{
			char const c = *a++;
			char const d = *--b;

			assert ( c == (d^3) );
		}
	}

	typedef libmaus2::autoarray::AutoArray < SortEntry > sort_array;
	typedef sort_array::unique_ptr_type sort_array_ptr;

	libmaus2::autoarray::AutoArray < sort_array_ptr > Asortarray(numthreads);
	for ( uint64_t i = 0; i < numthreads; ++i )
	{
		sort_array_ptr ptr(new sort_array);
		Asortarray[i] = std::move(ptr);
	}

	libmaus2::dazzler::align::AlignmentWriterArray AWA(tmpprefix + "_AWA",numthreads,tspace);

	#if defined(_OPENMP)
	#pragma omp parallel for schedule(dynamic,1) num_threads(numthreads)
	#endif
	for ( uint64_t i = 0; i < V.size(); ++i )
	{
		try
		{
			uint64_t const tid =
				#if defined(_OPENMP)
				omp_get_thread_num()
				#else
				0
				#endif
				;

			// libmaus2::lcs::NNPCor nnpcor;
			libmaus2::dazzler::align::AlignmentWriter & AW = AWA[tid];
			libmaus2::lcs::NNP nnpcor;
			libmaus2::lcs::NNPTraceContainer nnptrace;
			sort_array & Asort = *(Asortarray[tid]);
			uint64_t Asorto = 0;
			Coord const & repcoord = V[i];

			libmaus2::autoarray::AutoArray<libmaus2::lcs::NNPTraceContainer::DiagStrip> D;
			libmaus2::lcs::AlignmentTraceContainer ATC;

			uint64_t const indexseq = repcoord.rc ? (Pmeta->S.size()-repcoord.seq-1) : repcoord.seq;
			uint64_t const seqof = Pmeta->L[indexseq];
			uint64_t const seqlen = Pmeta->L[indexseq+1] - Pmeta->L[indexseq];

			uint64_t const basepos = repcoord.rc ? (seqlen - (repcoord.left + repcoord.length) - 1) : repcoord.left;

			char const * const pa = A.begin() + seqof + basepos;
			char const * const pe = pa + repcoord.length;

			{
			libmaus2::aio::StreamLock::lock_type::scope_lock_type slock(libmaus2::aio::StreamLock::cerrlock);
			std::cerr << repcoord << " " << indexseq << " " << seqof << std::endl;
			}

			uint64_t const l = pe-pa;
			uint64_t const numk = (l >= k) ? (l-k+1) : 0;

			// iterate over k-mers
			for ( uint64_t zzz = 0; zzz < numk; ++zzz )
			{
				char const * pc = pa + zzz;
				std::pair<uint64_t,uint64_t> const P = kcache.search(pc,k);

				// more than one hit?
				if ( P.second - P.first > 1 )
				{
					// iterate over hits
					for ( uint64_t zz = P.first; zz < P.second; ++zz )
					{
						// map coordinates
						uint64_t const p = GP.getPosition(zz);
						std::pair<uint64_t,uint64_t> const coL = cocache[p];
						std::pair<uint64_t,uint64_t> const coR = cocache[(p+k-1)%n];

						#if 0
						for ( uint64_t i = 0; i < k; ++i )
							std::cerr << libmaus2::fastx::remapChar(A[p+i]);
						std::cerr << "\t";
						for ( uint64_t i = 0; i < k; ++i )
							std::cerr << libmaus2::fastx::remapChar(pc[i]);
						std::cerr << std::endl;
						#endif

						if (
							// do not consider seeds spanning over a refseq boundary
							(coL.first == coR.first)
							&&
							// no self matches
							(A.begin() + p != pc)
						)
						{
							// push hit
							Asort.push(
								Asorto,
								SortEntry(
									basepos + zzz,
									coL.first,coL.second,k
								)
							);
						}
					}
				}
			}

			// sort hits
			std::sort(Asort.begin(),Asort.begin() + Asorto);

			// merge adjacent hits
			uint64_t ilow = 0;
			uint64_t o = 0;
			while ( ilow < Asorto )
			{
				uint64_t ihigh = ilow+1;
				while ( ihigh < Asorto && Asort[ilow].seq == Asort[ihigh].seq )
					++ihigh;

				// merge adjacent hits
				uint64_t dlow = ilow;
				while ( dlow < ihigh )
				{
					uint64_t dhigh = dlow+1;
					while ( dhigh < ihigh && Asort[dlow].getDiagonal() == Asort[dhigh].getDiagonal() )
						++dhigh;

					uint64_t alow = dlow;
					while ( alow < dhigh )
					{
						uint64_t ahigh = alow+1;

						while ( ahigh < dhigh && Asort[ahigh].overlap(Asort[ahigh-1]) )
							++ahigh;

						uint64_t const clow = Asort[alow].zzz;
						uint64_t const chigh = Asort[ahigh-1].zzz + Asort[ahigh-1].l;
						uint64_t const clen = chigh-clow;

						Asort[alow].l = clen;
						Asort[o++] = Asort[alow];

						// std::cerr << Asort[o] << std::endl;

						alow = ahigh;
					}


					dlow = dhigh;
				}

				ilow = ihigh;
			}
			Asorto = o;

			ilow = 0;

			uint64_t gen = 0;

			// iterate over hits
			while ( ilow < Asorto )
			{
				// same seq
				uint64_t ihigh = ilow+1;
				while ( ihigh < Asorto && Asort[ilow].seq == Asort[ihigh].seq )
					++ihigh;

				// split into hits sufficienlty close to each other
				std::vector < std::pair<uint64_t,uint64_t> > const R = gapFilter(Asort.begin(),ilow,ihigh,maxgap);

				// iterate over intervals
				for ( uint64_t i = 0; i < R.size(); ++i )
				{
					// interval high and low
					uint64_t glow = R[i].first;
					uint64_t ghigh = R[i].second;

					// sort by diagonal
					std::sort(
						Asort.begin()+glow,
						Asort.begin()+ghigh
					);

					uint64_t todo = ghigh-glow;

					// std::cerr << std::string(80,'*') << std::endl;
					while ( todo )
					{
						SortEntry & SE = Asort[glow];

						// std::cerr << SE << std::endl;

						#if 0
						char const * p0 = A.begin() + seqof + repcoord.left + SE.zzz;
						char const * p1 = A.begin() + Pmeta->L[SE.seq] + SE.pos;
						#endif

						// query sequence
						char const * p0a = A.begin() + Pmeta->L[indexseq+0];
						char const * p0e = A.begin() + Pmeta->L[indexseq+1];

						// ref sequence
						char const * p1a = A.begin() + Pmeta->L[SE.seq+0];
						char const * p1e = A.begin() + Pmeta->L[SE.seq+1];

						#if 1
						for ( uint64_t q = 0; q < SE.l; ++q )
						{
							assert ( p0a[SE.zzz + q] == p1a[SE.pos + q] );
						}
						#endif

						libmaus2::lcs::NNPAlignResult nnpres = nnpcor.align(
							p0a,p0e,SE.zzz,
							p1a,p1e,SE.pos,
							nnptrace,
							true /* check for self alignments */
						);

						// mark seeds touching alignment as inactive
						uint64_t const Do = nnptrace.getDiagStrips(nnpres.abpos,nnpres.bbpos,D);
						for ( uint64_t j = 0; j < Do; ++j )
						{
							libmaus2::lcs::NNPTraceContainer::DiagStrip const & DS = D[j];

							std::pair < SortEntry *, SortEntry * > const P =
								std::equal_range(
									Asort.begin()+glow,
									Asort.begin()+ghigh,
									SortEntry(DS.d, 0, 0, 0),
									SortEntryDiagComparator()
								);

							for ( SortEntry * p = P.first; p != P.second; ++p )
							{
								libmaus2::math::IntegerInterval<int64_t> const IA = p->getAntiDiagonalInterval();
								libmaus2::math::IntegerInterval<int64_t> const IB(
									DS.l,
									DS.h
								);

								assert ( DS.d == p->getDiagonal() );

								if ( !IA.intersection(IB).isEmpty() )
									p->active = false;
							}
						}


						// filter out inactive seeds
						SE.active = false;

						uint64_t o = glow;
						for ( uint64_t jj = glow; jj < ghigh; ++jj )
							if ( Asort[jj].active )
								Asort[o++] = Asort[jj];
						ghigh = o;
						todo = ghigh - glow;


						nnptrace.computeTrace(ATC);
						bool const alok = libmaus2::lcs::AlignmentTraceContainer::checkAlignment(
							ATC.ta,
							ATC.te,
							p0a + nnpres.abpos,
							p1a + nnpres.bbpos
						);
						assert ( alok );

						extend(p0a,p0e-p0a,p1a,p1e-p1a,nnpres,ATC,128 /* ex */,0.2);

						#if 0
						std::cerr << "[A] " << nnpres << std::endl;

						extendLeft(p0a,p1a,nnpres,ATC,128 /* ex */);
						std::cerr << "[B] " << nnpres << std::endl;

						extendRight(p0a,p0e-p0a,p1a,p1e-p1a,nnpres,ATC,128 /* ex */);
						std::cerr << "[C] " << nnpres << std::endl;
						#endif

						if ( nnpres.aepos - nnpres.abpos >= minreport )
						{
							//std::cerr << nnpres << " " << nnpres.aepos - nnpres.abpos << std::endl;

							uint64_t refseq = SE.seq;
							uint64_t qseq   = indexseq;

							// std::cerr << "refseq=" << refseq << " qseq=" << qseq << " numseq=" << numseq << std::endl;
							uint64_t const refseqlen = Pmeta->L[refseq+1] - Pmeta->L[refseq];
							uint64_t const qseqlen   = Pmeta->L[qseq+1  ] - Pmeta->L[qseq  ];

							if ( refseq >= numseq/2 )
							{
								//std::cerr << "inverting" << std::endl;


								#if 0
								assert ( nnpres.abpos <= nnpres.aepos );
								assert ( nnpres.aepos <= qseqlen );
								#endif

								std::swap(nnpres.abpos,nnpres.aepos);
								nnpres.abpos = qseqlen - nnpres.abpos;
								nnpres.aepos = qseqlen - nnpres.aepos;

								std::swap(nnpres.bbpos,nnpres.bepos);
								nnpres.bbpos = refseqlen - nnpres.bbpos;
								nnpres.bepos = refseqlen - nnpres.bepos;

								refseq = numseq - refseq - 1;
								qseq = numseq - qseq - 1;

								std::reverse(ATC.ta,ATC.te);
							}

							assert ( refseq < numseq/2 );

							ATC.swapRoles();

							uint64_t refpos = nnpres.bbpos;
							uint64_t reflen = nnpres.bepos-nnpres.bbpos;
							uint64_t qpos = nnpres.abpos;
							uint64_t qlen = nnpres.aepos-nnpres.abpos;

							char const * cref = A.begin() + Pmeta->L[refseq];
							char const * cq   = A.begin() + Pmeta->L[qseq  ];

							bool const alok = libmaus2::lcs::AlignmentTraceContainer::checkAlignment(
								ATC.ta,
								ATC.te,
								cref + refpos,
								cq + qpos
							);

							if ( ! alok )
							{
								libmaus2::lcs::AlignmentPrint::printAlignmentLines(
									std::cerr,
									cref + refpos,
									reflen,
									cq + qpos,
									qlen,
									80,
									ATC.ta,
									ATC.te,
									libmaus2::fastx::remapChar
								);
							}

							assert ( alok );

							bool const ovlrc = (qseq >= (numseq/2));
							uint64_t const ovlqseq = ovlrc ? (numseq - qseq - 1) : qseq;
							assert ( ovlqseq < numseq/2 );

							uint64_t ovlqpos_a = qpos       ;
							uint64_t ovlqpos_e = qpos + qlen;

							assert ( ovlqpos_a <= qseqlen );
							assert ( ovlqpos_e <= qseqlen );
							assert ( ovlqpos_a <= ovlqpos_e );

							#if 0
							if ( ovlrc )
							{
								std::swap(ovlqpos_a,ovlqpos_e);
								ovlqpos_a = qseqlen - ovlqpos_a;
								ovlqpos_e = qseqlen - ovlqpos_e;
							}
							#endif

							assert ( ovlqpos_a <= qseqlen );
							assert ( ovlqpos_e <= qseqlen );
							assert ( ovlqpos_a <= ovlqpos_e );

							#if 0
							{
								std::string sa(A.begin() + Pmeta->L[refseq]  + refpos   ,reflen);

								std::string rb(A.begin() + Pmeta->L[ovlqseq], Pmeta->L[ovlqseq+1]-Pmeta->L[ovlqseq]);
								if ( ovlrc )
									rb = libmaus2::fastx::reverseComplement(rb);
								std::string sb = rb.substr(ovlqpos_a,ovlqpos_e-ovlqpos_a);

								bool const alok = libmaus2::lcs::AlignmentTraceContainer::checkAlignment(
									ATC.ta,
									ATC.te,
									sa.begin(),
									sb.begin()
								);

								if ( ! alok )
								{
									std::cerr << ovlrc << std::endl;
								}
								assert(alok);
							}
							#endif

							libmaus2::dazzler::align::Overlap const OVL = libmaus2::dazzler::align::Overlap::computeOverlap(
								ovlrc,
								refseq,
								ovlqseq,
								refpos,
								refpos + reflen,
								ovlqpos_a,
								ovlqpos_e,
								tspace,
								ATC
							);

							#if 0
							std::cerr << OVL << std::endl;

							{
								std::string sa(
									A.begin() + Pmeta->L[OVL.aread],
									Pmeta->L[OVL.aread+1]-Pmeta->L[OVL.aread]
								);
								std::string sb(
									A.begin() + Pmeta->L[OVL.bread],
									Pmeta->L[OVL.bread+1]-Pmeta->L[OVL.bread]
								);
								if ( OVL.isInverse() )
									sb = libmaus2::fastx::reverseComplement(sb);

								libmaus2::lcs::NP np;
								libmaus2::lcs::AlignmentTraceContainer ATC;
								OVL.computeTrace(
									reinterpret_cast<uint8_t const *>(sa.c_str()),
									reinterpret_cast<uint8_t const *>(sb.c_str()),
									tspace,
									ATC,
									np
								);

								std::cerr << ATC.getAlignmentStatistics() << std::endl;
							}
							#endif

							//std::cerr << OVL << std::endl;

							AW.put(OVL);

							if ( (++gen % 1024) == 0 )
							{
								libmaus2::aio::StreamLock::lock_type::scope_lock_type slock(libmaus2::aio::StreamLock::cerrlock);
								std::cerr << "[tid=" << tid << "] " << V[i] << " seq=" << Asort[ilow].seq << " proc=" << i << "/" << R.size() << std::endl;
							}

						}

					}
				}

				ilow = ihigh;
			}
		}
		catch(std::exception const & ex)
		{
			std::cerr << ex.what() << std::endl;
			std::terminate();
		}
	}

	AWA.merge(outfn,tmpprefix+"_tmp");

	return EXIT_SUCCESS;
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
		else if ( arg.uniqueArgPresent("h") || arg.uniqueArgPresent("help") || arg.size() < 2 )
		{
			std::cerr << "This is " << PACKAGE_NAME << " version " << PACKAGE_VERSION << "." << std::endl;
			std::cerr << PACKAGE_NAME << " is distributed under version 3 of the GPL." << std::endl;
			std::cerr << "\n";
			std::cerr << "usage: " << arg.progname << " [options] out.las ref.fasta ref.fasta.selfie\n";
			std::cerr << "\n";
			std::cerr << "The following options can be used (no space between option name and parameter allowed):\n\n";
			//std::cerr << helpMessage(arg);
			return EXIT_SUCCESS;
		}
		else
		{
			return remapselfie(arg);
		}
	}
	catch(std::exception const & ex)
	{
		std::cerr << ex.what() << std::endl;
		return EXIT_FAILURE;
	}

}
