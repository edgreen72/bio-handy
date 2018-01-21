# bio-handy
This repo contains random scripts for doing things that are occassionally useful in sequence and genome analysis.

## MiSeqRunCheck.pl
```
MiSeqRunCheck.pl -s <SampleSheet.csv> -d <DemultiplexSummaryF1L1.txt
Reports the IDs and counts of the most common indeces seen in a run
and writes NA if the index was not listed in the SampleSheet.csv
```

## align2consensus.pl
```
Takes as input a multiple-sequence alignment file in clustalw
format. Generates a single, fasta consensus sequence by
calling the most common base at each postion. Where there
is no consensus over the user-defined percent identity, the
output sequence contains an N.
Requires Bio::AlignIO and Bio::SeqIO from bioperl
  -a <alignment file; clustalw format>
  -c <percent ID for consensus call; default 50>
  -I <Identifier to use; default Consensus>
```

## bam-bc-count.pl
```
bam-bc-count.pl -b <bam file> -5 <P5 barcode> -7 <P7 barcodes>
                -q <map quality cutoff>
Reports a count of the observed number of each barcode combination
in the input bam file. Barcodes are read from the BC field.
Requires BIo::DB::Sam
```

## count-restriction-sites

## fasta-gc-window.pl
```
fasta-gc-window.pl -f <fasta file> -w <window; DEF = 100000>
Makes a table of:
1. Sequence ID
2. Window start position
3. Percent GC in window
Requires Bio::SeqIO
```

## fastq-barcode-split.pl
```
fastq-barcode-split.pl -f <forward fastq file> -r <reverse fastq file> -l [BARCODE length]
                       -o <root name for output file(s) to make>
The program reads one or a pair of fastq files specified by the -f and -r options.
If no -r file is given, it's ass-u-med that the single input file has reads with
no adapters and barcodes at the very beginning and end of the reads.
If both -f and -r files are given, they are ass-u-med to be the forward and reverse
reads from a paired-end run. Reads must be in the same order.
It takes the first -l bases of the forward and reverse reads and treats them as a
composite barcode. It reports the number of read pairs with each observed barcode.
For example, a read pair like this, using -l 5
@M00160:20:000000000-ANYJ6:1:1101:10222:1120 1:N:0:8\nGATGAGGCTTAGAAGCAAGATCGGAAGAGCACACGTCTGAACTCCAGTCACTTGGTCAACATCGTTTGCCGTCTTC
+
CC@CCGGGEFFACFFFGC6,,CCCCC+6CE@<FFDCFFG9C@F9<@CFDAFCG<,;CAF,,;@,C,C,EF+CEF7C\n[other read]
@M00160:20:000000000-ANYJ6:1:1101:10222:1120 2:N:0:8\nTGCTTCTAAGCCTAATCAGATCGGAAGAGCGTCGTGTAGNGAAAGANTNNAGNTNTCGNTNGGNGCCGTATCATTA
+
CCCCCGGGGG?8F,;CFFEDFF,@FGD8@D@F,66BCEC#6,CC@<#:##,,#6#::C#,#:,#:C,@C@,,<,,<
Would be classified as GATGA.TGCTT
```

## orf-scan.pl
```
orf-scan.pl Version 1 -f <fasta genome> -s <STARTs> -t <TERMs> -m <MIN>
Scans the fasta genome sequence for the input start and termination
codons and outputs all observed ORFs greater than MIN length.
STARTs and TERMs should be a colon delimited list. For example:
-s ATG -t TAG:TAA:TGA
Writes an output file with the folowing columns
1. ID
2. Strand (Forward or Reverse)
3. Frame (0, 1, or 2)
4. Position of first base of start codon
5. Position of first base of stop codon
-s DEFAULT = ATG
-t DEFAULT = TAG:TAA:TGA
-m DEFAULT = 27
```

## longest-orf.pl
```
Takes input from STDIN, assumed to output of orf-scan.pl and
finds the longest-orf.pl
```

## pss-bam.pl
```
pss-bam.pl v 0.06 -f <fasta file> -b <bam file>
           -r <region length; default = 15
           -l <min length filter; default = 0>
           -L <max length filter; default = 250000000>
           -q <map quality filter>
           -F <front barcode sequence>
           -B <back barcode sequence>
           -m <run in merged mode => reads are merged>
Uses the Bio::DB::Sam module to make map-damage like data
suitable for plotting in gnuplot.
```