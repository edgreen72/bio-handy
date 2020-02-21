#CC=gcc
CFLAGS=-O2
#CFLAGS=-gdwarf-2 -g

fasta-genome-io.o : fasta-genome-io.h fasta-genome-io.c
	echo "Making fasta-genome-io.o..."
	$(CC) $(CFLAGS) fasta-genome-io.c -c -lz -o fasta-genome-io.o

test-fasta-genome : test-fasta-genome.c fasta-genome-io.o
	echo "Making test-fasta-genome..."
	$(CC) $(CFLAGS) fasta-genome-io.o test-fasta-genome.c -lz -o test-fasta-genome

fastq-io.o : fastq-io.h fastq-io.c
	echo "Making fastq-io.o ..."
	$(CC) $(CFLAGS) fastq-io.c -c -lz -o fastq-io.o

fastq-dinuc-count : fastq-dinuc-count.c fastq-io.o
	echo "Making fastq-dinuc-count..."
	$(CC) $(CFLAGS) fastq-io.o fastq-dinuc-count.c -lz -o fastq-dinuc-count 

what-adapter : what-adapter.c fastq-io.o
	echo "Making what-adapter..."
	$(CC) $(CFLAGS) fastq-io.o what-adapter.c -lz -o what-adapter
