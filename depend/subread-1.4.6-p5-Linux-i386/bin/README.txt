Installation
--------------
You may use Subread package by directly downloading a binary release suitable for your operating system (no compilation is needed), or you may build it from the source. Here we describe how to install it from the source.

Download the latest version of Subread package from http://subread.sourceforge.net/. The source release includes a keyword 'source' in the name of the tar ball. Uncompress the tar ball, enter the 'src' directory and issue the following command to build it for Linux OS :

make -f Makefile.Linux

For Mac OS, use command:

make -f Makefile.MacOS

For FreeBSD OS, use command:

gmake -f Makefile.FreeBSD

For OpenSolaris/Oracle Solaris OS, use command:

gmake -f Makefile.SunOS

If the build is successful, a new directory called 'bin' will be created under the home directory of the package (ie. one level up from 'src' directory). The 'bin directory contains all the generated executables. To enable easy access to these executables, you may copy the executables to a system directory such as '/usr/bin' or add the path to the executables to your search path (add path to your environment variable `PATH').

Content
--------------
annotation	Directory including NCBI RefSeq gene annotations for genomes 'hg19', 'mm10' and 'mm9'.
			Each row is an exon. Entrez gene identifiers and chromosomal coordinates are provided for each exon.
bin			Directory including executables after compilation (or directly available from a binary release). 
doc			Directory including the users manual.
LICENSE		The license agreement for using this package.
README.txt	This file.
src			Directory including source code (binary releases do not have this directory).
test		Directory including test data and scripts.

A Quick Start
--------------
An index should be built before carrying out read alignments: 

subread-buildindex -o my_index chr1.fa chr2.fa ...
(You may provide a single FASTA file including all chromosomal sequences).

With built index, you can now align reads to the reference genome. Align single-end reads: 

subread-align -i my_index -r reads.txt -o subread_results.sam

Align paired-end reads:

subread-align -i my_index -r reads1.txt -R reads2.txt -o subread_results_PE.sam

Detect exon-exon junctions from RNA-seq data (read mapping results are also generated):

subjunc -i my_index -r reads1.txt -R reads2.txt -o subjunc_results.sam

Assign mapped reads to genomic features (eg. genes):

featureCounts -a annotation.gtf -o counts.txt subread_results.sam

Tutorials
-------------------
A short tutorial for Subread - http://bioinf.wehi.edu.au/subread
A short tutorial for Subjunc - http://bioinf.wehi.edu.au/subjunc
A short tutorial for featureCounts - http://bioinf.wehi.edu.au/featureCounts

Users Guide
--------------
Users Guide can be found in the 'doc' subdirectory. It provides comprehensive descriptions to the programs included in this package.

Citation
--------------
If you use Subread or Subjunc aligners, please cite:

Liao Y, Smyth GK and Shi W. The Subread aligner: fast, accurate and scalable read mapping by seed-and-vote. Nucleic Acids Research, 41(10):e108, 2013

If you use the featureCounts program, please cite:

Liao Y, Smyth GK and Shi W. featureCounts: an efficient general-purpose program for assigning sequence reads to genomic features. Bioinformatics, 2013. doi: 10.1093/bioinformatics/btt656

Get help
--------------
You may subscribe to the SeqAnswer forum (http://www.seqanswers.com) or the Bioconductor mailing list (http://bioconductor.org/) to get help. Alternatively, you may directly contact Wei Shi (shi at wehi dot edu dot au) or Yang Liao (liao at wehi dot edu dot au) for help.
