# yalora
Yalora is a prototypical long read aligner. The name is an acronym for **Yet another long read aligner**.
It does not have the speed of other aligners particularly on larger repetitive genomes but is constructed to be reasonably accurate/sensitive.

## Source

The yalora source code is hosted on github:

	git@github.com:gt1/yalora.git

## Compilation of yalora

yalora needs libmaus2 [https://github.com/gt1/libmaus2] . When libmaus2
is installed in ${LIBMAUSPREFIX} then yalora can be compiled and
installed in ${HOME}/yalora using

	- autoreconf -i -f
	- ./configure --with-libmaus2=${LIBMAUSPREFIX} \
		--prefix=${HOME}/yalora
	- make install

## Running alignments

The aligner can be run using

```
yalora [options] ref.fasta < reads.fasta > reads_aligned.bam
```

with reference and reads in the FastA [https://en.wikipedia.org/wiki/FASTA] format. It outputs alignments in the bam [https://samtools.github.io/hts-specs/SAMv1.pdf] format, by default.
The output format can be set to sam or cram instead (where cram is most likely not a good idea as the output in general is not sorted by coordinate but in input order).

Yalora may output more than one alignment per read if it finds several alignments of similar quality. It will mark one of the alignments with maximal alignment score as primary and all other alignments
for a read as secondary.

Tests with reads simulated from real genomes (E.coli, D.mel and HG19) at 10% error rates show that the primary alignment is usually the one mapped in the correct place, however for reads originating
from repetitive regions the correct alignment may not be the one marked as primary. In most cases the correct region to place a read is contained in the aligners output though.

The seed selection strategy shares some similarity with BWA's [https://arxiv.org/abs/1303.3997] MEM algorithm, the following stages are different though.

## Options

Yalora has the following options (no space between option and argument allowed):

* -i: input format, e.g. `-ifasta`, only fasta supported so far
* -o: output format, e.g. `-obam`, bam or sam
* -t: number of threads, e.g. `-t8`, by default this is chosen as the number of logical CPU cores on the executing machine
* --constructmem: memory guideline for constructing the index used (BWT+sampled suffix array) in bytes, e.g. `--constructmem16g`, by default this is chosen is 3/4 of the physical RAM detected on the executing machine
* --sasamplingrate: suffix array sampling rate, e.g. `--sasamplingrate16`, 32 by default
* -T: prefix for temporary files produced during index construction and output reordering, e.g. `-T/path/to/dir/tmp_`, by default will generate files in the current directory with names based on the executable name, executing machine's host name, process ID and timestamp
* --chainminscore: minimum score required for seed chains, e.g. `--chainminscore20`, default 20. Chains with smaller scores are discarded. Increasing this values will increase speed and decrease sensitivity.
* --minlen: minimum seed length, e.g. `--minlen20`, default 20. Matches shorter than this length are ignored.
* --minsplitlength: minimum seed length to consider splitting long seeds, e.g. `--minsplitlength28`, default 28. This option is very similar to the reseeding described in section 2.1.1 of the BWA MEM paper [https://arxiv.org/abs/1303.3997].
* --minsplitsize: maximum frequency to consider splitting, e.g. `--minsplitsize10`, default 10.
* -K: kmer cache K size, e.g. `-K12`, default 12.
* --maxocc: maximum frequency of seeds used, e.g. `--maxocc500`, default 500. Seeds more frequent then this threshold are subsampled.
* --minfreq: minimum frequency required for seeds, e.g. `--minfreq1`, default 1.
* --limit: maximum seed extension left and right, e.g. `--limit32`, default 32. Maximum number of bases extended to left and right starting from a central position when looking for seeds.
* --maxxdist: maximum distance allowed for chaining seeds in any of the two coordinates, e.g. `--maxxdist1000`, default 1000.
* --activemax: number of chains required to cover/dominate another chain until it is discarded, e.g. `--activemax1`, default 1.
* --fracmul: fragment overlap threshold counter, e.g. `--fracmul95`, default 95. fracmul/fracdiv designates the fraction of bases two chains need to overlap by to consider them for testing whether one dominates the other
* --fracdiv: fragment overlap threshold denominator, e.g. `--fracdiv100`, default 100.
* --chaindommul: chain domination threshold counter, e.g. `--chaindommul95`, default 95. A chain covers/dominates another chain if it overlaps fracmul/fracdiv of its bases and its score multiplied by chaindommul/chaindomdiv exceeds the other chains score/range.
* --chaindomdiv: chain domination threshold denominator, e.g. `--chaindomdiv100`, default 100.
* --algndommul: alignment domination threshold counter, e.g. `--algndommul95`, default 95. An alignment covers/dominates another alignment if it overlaps fracmul/fracdiv of its bases and its score multiplied by algndommul/algndomdiv exceeds the other alignments score/range.
* --algndomdiv: alignment domination threshold denominator, e.g. `--algndomdiv100`, default 100.

The default values should be reasonable for aligning 3rd generation long read data.
