/*
    libmaus2
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
#include <libmaus2/lcs/SMEMProcessor.hpp>

#include <libmaus2/fastx/DNAIndexMetaDataBigBandBiDir.hpp>
#include <libmaus2/fastx/FastAToCompact4BigBandBiDir.hpp>
#include <libmaus2/util/ArgInfo.hpp>
#include <libmaus2/util/ArgParser.hpp>
#include <libmaus2/suffixsort/bwtb3m/BwtMergeSortOptions.hpp>
#include <libmaus2/suffixsort/bwtb3m/BwtMergeSort.hpp>
#include <libmaus2/aio/DebugLineOutputStream.hpp>
#include <libmaus2/sorting/SortingBufferedOutputFile.hpp>

static std::vector < std::string > stateVec;
static std::vector < std::string > printVec;
static libmaus2::parallel::StdSpinLock stateVecLock;

static void setState(uint64_t const tid, std::string const & s)
{
	libmaus2::parallel::ScopeStdSpinLock slock(stateVecLock);
	stateVec[tid] = s;
}

static void printState()
{
	libmaus2::aio::DebugLineOutputStream DLOS(std::cerr,libmaus2::aio::StreamLock::cerrlock);
	libmaus2::parallel::ScopeStdSpinLock slock(stateVecLock);

	DLOS << static_cast<char>(27) << 'c';
	for ( uint64_t i = 0; i < stateVec.size(); ++i )
	{
		DLOS << i << "\t" << stateVec[i] << "\n";
	}
}


struct CoordinatePair
{
	libmaus2::fastx::DNAIndexMetaDataBigBandBiDir::Coordinates A;
	libmaus2::fastx::DNAIndexMetaDataBigBandBiDir::Coordinates B;
	libmaus2::lcs::NNPAlignResult res;

	CoordinatePair() {}
	CoordinatePair(
		libmaus2::fastx::DNAIndexMetaDataBigBandBiDir::Coordinates const & rA,
		libmaus2::fastx::DNAIndexMetaDataBigBandBiDir::Coordinates const & rB,
		libmaus2::lcs::NNPAlignResult const & rres
	) : A(rA), B(rB), res(rres)
	{

	}
	CoordinatePair(std::istream & in)
	{
		deserialise(in);
	}

	void deserialise(std::istream & in)
	{
		A.deserialise(in);
		B.deserialise(in);
		res.deserialise(in);
	}

	void serialise(std::ostream & out) const
	{
		A.serialise(out);
		B.serialise(out);
		res.serialise(out);
	}

	bool operator<(CoordinatePair const & C) const
	{
		if ( A < C.A )
			return true;
		else if ( C.A < A )
			return false;
		else if ( B < C.B )
			return true;
		else if ( C.B < B )
			return false;
		else
			return false;
	}
};

struct GammaInterval
{
	uint64_t first;
	uint64_t second;

	GammaInterval() {}
	GammaInterval(uint64_t const rfirst, uint64_t const rsecond) : first(rfirst), second(rsecond) {}
	GammaInterval(std::istream & in)
	{
		deserialise(in);
	}

	void deserialise(std::istream & in)
	{
		first = libmaus2::util::NumberSerialisation::deserialiseNumber(in);
		second = libmaus2::util::NumberSerialisation::deserialiseNumber(in);
	}

	void serialise(std::ostream & out) const
	{
		libmaus2::util::NumberSerialisation::serialiseNumber(out,first);
		libmaus2::util::NumberSerialisation::serialiseNumber(out,second);
	}

	int64_t getDiam() const
	{
		return static_cast<int64_t>(second)-static_cast<int64_t>(first);
	}

	bool operator<(GammaInterval const & O) const
	{
		if ( getDiam() != O.getDiam() )
			return getDiam() > O.getDiam();
		else if ( first != O.first )
			return first < O.first;
		else
			return second < O.second;
	}
};

void selfie(libmaus2::util::ArgParser const & arg, std::string const & fn)
{
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

	uint64_t const n = Prank->size();
	libmaus2::autoarray::AutoArray<char> A(n,false);
	libmaus2::bitio::CompactDecoderWrapper CDW(compactfn);
	CDW.read(A.begin(),n);
	assert ( CDW.gcount() == static_cast<int64_t>(n) );

	uint64_t const minfreq = 2;
	uint64_t const minlen = 20;
	uint64_t const limit = 32;
	uint64_t const minsplitlength = 28;
	uint64_t const minsplitsize = 10;
	uint64_t const maxxdist = 1000;
	uint64_t const activemax = 1;
	uint64_t const fracmul = 95;
	uint64_t const fracdiv = 100;
	bool const selfcheck = true;
	uint64_t const chainminscore = arg.uniqueArgPresent("chainminscore") ? arg.getUnsignedNumericArg<uint64_t>("chainminscore") : 20;
	uint64_t const maxocc = 500;
	uint64_t const minprintlength = 1024;
	uint64_t const algndommul = 95;
	uint64_t const algndomdiv = 100;
	uint64_t const chaindommul = 95;
	uint64_t const chaindomdiv = 100;
	double const maxerr = arg.uniqueArgPresent("maxerr") ? arg.getParsedArg<double>("maxerr") : std::numeric_limits<double>::max();

	uint64_t const cachek = arg.uniqueArgPresent("K") ? arg.getUnsignedNumericArg<uint64_t>("K") : 12;
	uint64_t const maxpacksize = arg.uniqueArgPresent("P") ? arg.getUnsignedNumericArg<uint64_t>("P") : 128ull*1024ull*1024ull;
	std::cerr << "[V] generating " << cachek << "-mer cache...";
	libmaus2::rank::DNARankKmerCache::unique_ptr_type Pcache(new libmaus2::rank::DNARankKmerCache(*Prank,cachek,numthreads));
	std::cerr << "done." << std::endl;

	std::string const deftmp = libmaus2::util::ArgInfo::getDefaultTmpFileName(arg.progname);
	libmaus2::util::TempFileNameGenerator tmpgen(deftmp,3);

	std::string const sorttmp = tmpgen.getFileName();
	libmaus2::util::TempFileRemovalContainer::addTempFile(sorttmp);
	libmaus2::sorting::SortingBufferedOutputFile<CoordinatePair> CPS(sorttmp);
	libmaus2::parallel::StdSpinLock CPSlock;

	uint64_t acc_s = 0;
	for ( uint64_t zz = 0; zz < Pmeta->S.size(); )
	{
		uint64_t zze = zz;
		uint64_t pack_s = Pmeta->S[zze++].l;

		while ( zze < Pmeta->S.size() && pack_s + Pmeta->S[zze].l <= maxpacksize )
			pack_s += Pmeta->S[zze++].l;

		// std::cerr << "[V] " << zz << "-" << zze << " pack_s=" << pack_s << std::endl;

		zz = zze;

		uint64_t const low = acc_s;
		uint64_t const high = acc_s + pack_s;

		std::cerr << "[V] low=" << low << " high=" << high << " acc_s=" << acc_s << " pack_s=" << pack_s << std::endl;

		std::string const activefn =
			libmaus2::rank::DNARankSMEMComputation::activeParallel(tmpgen,*Pcache,A.begin(),low,high,minfreq,minlen,numthreads,maxxdist + 2*(minlen-1));

		libmaus2::gamma::GammaIntervalDecoder::unique_ptr_type Pdec(new libmaus2::gamma::GammaIntervalDecoder(std::vector<std::string>(1,activefn),0/*offset */,1 /* numthreads */));

		std::string const sortinfn = tmpgen.getFileName(true);
		libmaus2::sorting::SerialisingSortingBufferedOutputFile<GammaInterval>::unique_ptr_type sptr(
			new libmaus2::sorting::SerialisingSortingBufferedOutputFile<GammaInterval>(sortinfn)
		);

		{
			std::pair<uint64_t,uint64_t> P;
			while ( Pdec->getNext(P) )
			{
				sptr->put(
					GammaInterval(P.first,P.second)
				);
			}
		}

		libmaus2::sorting::SerialisingSortingBufferedOutputFile<GammaInterval>::merger_ptr_type Pmerger(
			sptr->getMerger()
		);

		struct LockedGet
		{
			libmaus2::parallel::StdSpinLock lock;
			// libmaus2::gamma::GammaIntervalDecoder& dec;
			libmaus2::sorting::SerialisingSortingBufferedOutputFile<GammaInterval>::merger_ptr_type & Pmerger;

			LockedGet(libmaus2::sorting::SerialisingSortingBufferedOutputFile<GammaInterval>::merger_ptr_type & rPmerger) : Pmerger(rPmerger)
			{
			}

			bool getNext(std::pair<uint64_t,uint64_t> & P)
			{
				bool ok = false;
				{
					libmaus2::parallel::ScopeStdSpinLock slock(lock);
					GammaInterval Q;
					ok = Pmerger->getNext(Q);
					if ( ok )
					{
						P.first = Q.first;
						P.second = Q.second;
					}
				}
				return ok;
			}
		};

		libmaus2::autoarray::AutoArray < std::pair<uint64_t,uint64_t> > VP(numthreads);
		LockedGet LG(Pmerger);

		libmaus2::fastx::CoordinateCacheBiDir cocache(*Prank,*Pmeta,16 /* blockshfit */);

		typedef libmaus2::suffixsort::bwtb3m::BwtMergeSortResult::BareSimpleSampledSuffixArray sa_type;
		typedef libmaus2::lcs::SMEMProcessor<sa_type> smem_proc_type;
		libmaus2::autoarray::AutoArray < smem_proc_type::unique_ptr_type > Aproc(numthreads);
		for ( uint64_t i = 0; i < numthreads; ++i )
		{
			smem_proc_type::unique_ptr_type proc(new smem_proc_type(
				*Pmeta,cocache,*Prank,BSSSA,A.begin(),maxxdist,activemax,fracmul,fracdiv,selfcheck,chainminscore,maxocc,algndommul,algndomdiv,chaindommul,chaindomdiv,
				libmaus2::lcs::NNP::getDefaultMaxWindowError(),libmaus2::lcs::NNP::getDefaultMaxBack(),false /* domsameref */
				)
			);
			Aproc[i] = std::move(proc);
		}

		stateVec.resize(numthreads);
		for ( uint64_t i = 0; i < numthreads; ++i )
			setState(i,"idle");


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
			std::pair<uint64_t,uint64_t> & P = VP[tid];
			smem_proc_type & proc = *(Aproc[tid]);

			struct SelfieVerbosity : public smem_proc_type::Verbosity
			{
				uint64_t tid;
				std::string prefix;

				SelfieVerbosity(uint64_t const rtid, std::string const & rprefix)
				: tid(rtid), prefix(rprefix)
				{

				}

				void operator()(libmaus2::rank::DNARankMEM const & smem, uint64_t const z) const
				{
					std::ostringstream ostr;
					ostr << prefix << "\t" << z << "\t" << smem;
					setState(tid,ostr.str());
					printState();
				}
			};

			while ( LG.getNext(P) )
			{
				uint64_t const smemleft = std::max(static_cast<int64_t>(0),static_cast<int64_t>(P.first)-static_cast<int64_t>(minlen-1));
				uint64_t const smemright = std::min(P.second+minlen,n);

				std::ostringstream msgstr;
				msgstr << "[" << smemleft << "," << smemright << ")";

				setState(tid,msgstr.str());
				printState();

				libmaus2::rank::DNARankSMEMComputation::SMEMEnumerator<char const *> senum(
					*Prank,A.begin(),
					smemleft,
					smemright,
					minfreq,
					minlen,
					limit,
					minsplitlength,
					minsplitsize);

				SelfieVerbosity SV(tid,msgstr.str());

				proc.process(senum,A.begin(),n,minprintlength,maxerr,SV);
				// proc.printAlignments(minprintlength);

				std::pair<libmaus2::lcs::ChainAlignment const *, libmaus2::lcs::ChainAlignment const *> const AP =
					proc.getAlignments();

				for ( libmaus2::lcs::ChainAlignment const * it = AP.first; it != AP.second; ++it )
				{
					libmaus2::lcs::ChainAlignment const & CA = *it;
					libmaus2::lcs::NNPAlignResult const & res = CA.res;

					std::vector<libmaus2::fastx::DNAIndexMetaDataBigBandBiDir::Coordinates> const VA = Pmeta->mapCoordinatePairToList(res.abpos,res.aepos);
					std::vector<libmaus2::fastx::DNAIndexMetaDataBigBandBiDir::Coordinates> const VB = Pmeta->mapCoordinatePairToList(res.bbpos,res.bepos);

					if ( VA.size() == 1 && VB.size() == 1 )
					{
						CoordinatePair CP(VA[0],VB[0],res);
						libmaus2::parallel::ScopeStdSpinLock slock(CPSlock);
						CPS.put(CP);
					}
				}

				setState(tid,"idle");
				printState();

				#if 0
				std::cerr << "P=[" << P.first << "," << P.second << ")" << std::endl;

				{
					std::vector<libmaus2::rank::DNARankMEM> SMEM;
					libmaus2::rank::DNARankSMEMComputation::smemLimitedParallel(
						*Prank,
						*Pcache,
						A.begin(),
						P.first,
						P.second,
						n,
						minfreq,
						minlen,
						limit,
						SMEM,
						1 /* threads */);
					std::cerr << "[V] number of SMEMs is " << SMEM.size() << std::endl;

					// deallocate k-mer cache
					// Pcache.reset();

					std::vector<libmaus2::rank::DNARankMEM> SMEMsplit;
					libmaus2::rank::DNARankSMEMComputation::smemLimitedParallelSplit(*Prank,A.begin(),P.first,P.second,minlen,limit,minsplitlength,minsplitsize,SMEM,SMEMsplit,1 /* threads */);
					std::cerr << "[V] number of split SMEMs is " << SMEMsplit.size() << std::endl;

					// insert split SMEMs into regular SMEMs
					std::copy(SMEMsplit.begin(),SMEMsplit.end(),std::back_insert_iterator< std::vector<libmaus2::rank::DNARankMEM> >(SMEM));
					//libmaus2::sorting::InPlaceParallelSort::inplacesort2(SMEM.begin(),SMEM.end(),numthreads,libmaus2::rank::DNARankMEMPosComparator());
					std::sort(SMEM.begin(),SMEM.end(),libmaus2::rank::DNARankMEMPosComparator());

					SMEM.resize(std::unique(SMEM.begin(),SMEM.end())-SMEM.begin());

					libmaus2::rank::DNARankSMEMComputation::SMEMEnumerator<char const *> senum(
						*Prank,A.begin(),
						std::max(static_cast<int64_t>(0),static_cast<int64_t>(P.first)-static_cast<int64_t>(minlen-1)),
						std::min(P.second+minlen,n),
						minfreq,
						minlen,
						limit,
						minsplitlength,
						minsplitsize);

					libmaus2::rank::DNARankMEM smem;
					uint64_t c = 0;
					while ( senum.getNext(smem) )
					{
						// std::cerr << "ccc=" << smem << std::endl;

						if ( c >= SMEM.size() || smem != SMEM[c] )
						{
							std::cerr << "mismatch " << c << " " << smem;
							if ( c < SMEM.size() )
								std::cerr << " != " << SMEM[c];
							else
								std::cerr << " ???";
							std::cerr << std::endl;
						}
						else
						{
							std::cerr << "match " << c << " " << smem << " " << SMEM[c] << std::endl;
						}

						++c;
					}

					std::cerr << "c=" << c << " V=" << SMEM.size() << std::endl;
				}
				#endif
			}
		}

		acc_s += pack_s;
	}

	libmaus2::sorting::SortingBufferedOutputFile<CoordinatePair>::merger_ptr_type Pmerger(CPS.getMerger());
	CoordinatePair CP;
	while ( Pmerger->getNext(CP) )
	{
		std::ostringstream ostr;
		CP.A.print(ostr);
		ostr << " ";
		CP.B.print(ostr);
		ostr << " ";
		ostr << CP.res.getErrorRate();

		std::cout << ostr.str() << std::endl;
	}
}

int main(int argc, char * argv[])
{
	try
	{
		libmaus2::util::ArgParser arg(argc,argv);
		for ( uint64_t i = 0; i < arg.size(); ++i )
			selfie(arg,arg[i]);
	}
	catch(std::exception const & ex)
	{
		std::cerr << ex.what() << std::endl;
		return EXIT_FAILURE;
	}
}
