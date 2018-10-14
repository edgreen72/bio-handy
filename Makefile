CC=gcc
#CFLAGS=-gdwarf-2 -g
CFLAGS=-O2

fastq-io.o : fastq-io.h fastq-io.c
	echo "Making fastq-io.o ..."
	$(CC) $(CFLAGS) fastq-io.c -c -lz -o fastq-io.o

fastq-dinuc-count : fastq-dinuc-count.c fastq-io.o
	echo "Making fastq-dinuc-count..."
	$(CC) $(CFLAGS) fastq-io.o fastq-dinuc-count.c -lz -o fastq-dinuc-count 
