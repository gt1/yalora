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

static std::vector < std::string > stateVec;
static std::vector < std::string > printVec;
static libmaus2::parallel::PosixSpinLock stateVecLock;

static void setState(uint64_t const tid, std::string const & s)
{
	libmaus2::parallel::ScopePosixSpinLock slock(stateVecLock);
	stateVec[tid] = s;
}

static void printState()
{
	libmaus2::aio::DebugLineOutputStream DLOS(std::cerr,libmaus2::aio::StreamLock::cerrlock);
	libmaus2::parallel::ScopePosixSpinLock slock(stateVecLock);

	for ( uint64_t i = 0; i < stateVec.size(); ++i )
	{
		DLOS << i << "\t" << stateVec[i] << "\n";
	}
}

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

		std::string const activefn =
			libmaus2::rank::DNARankSMEMComputation::activeParallel(tmpgen,*Pcache,A.begin(),low,high,minfreq,minlen,numthreads,maxxdist + 2*(minlen-1));

		libmaus2::gamma::GammaIntervalDecoder::unique_ptr_type Pdec(new libmaus2::gamma::GammaIntervalDecoder(std::vector<std::string>(1,activefn),0/*offset */,1 /* numthreads */));

		struct LockedGet
		{
			libmaus2::parallel::PosixSpinLock lock;
			libmaus2::gamma::GammaIntervalDecoder& dec;

			LockedGet(libmaus2::gamma::GammaIntervalDecoder & rdec) : dec(rdec)
			{
			}

			bool getNext(std::pair<uint64_t,uint64_t> & P)
			{
				bool ok = false;
				{
					libmaus2::parallel::ScopePosixSpinLock slock(lock);
					ok = dec.getNext(P);
				}
				return ok;
			}
		};

		libmaus2::autoarray::AutoArray < std::pair<uint64_t,uint64_t> > VP(numthreads);
		LockedGet LG(*Pdec);

		libmaus2::fastx::CoordinateCacheBiDir cocache(*Prank,*Pmeta,16 /* blockshfit */);

		typedef libmaus2::suffixsort::bwtb3m::BwtMergeSortResult::BareSimpleSampledSuffixArray sa_type;
		typedef libmaus2::lcs::SMEMProcessor<sa_type> smem_proc_type;
		libmaus2::autoarray::AutoArray < smem_proc_type::unique_ptr_type > Aproc(numthreads);
		for ( uint64_t i = 0; i < numthreads; ++i )
		{
			smem_proc_type::unique_ptr_type proc(new smem_proc_type(
				*Pmeta,cocache,*Prank,BSSSA,A.begin(),maxxdist,activemax,fracmul,fracdiv,selfcheck,chainminscore,maxocc,algndommul,algndomdiv,chaindommul,chaindomdiv,
				libmaus2::lcs::NNP::getDefaultMaxWindowError(),libmaus2::lcs::NNP::getDefaultMaxBack()
				)
			);
			Aproc[i] = UNIQUE_PTR_MOVE(proc);
		}

		stateVec.resize(numthreads);

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

				proc.process(senum,A.begin(),n,minprintlength,maxerr);
				proc.printAlignments(minprintlength);

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
