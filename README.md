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
