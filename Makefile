CC=gcc
CFLAGS=-O2
#CFLAGS=-gdwarf-2 -g


fastq-io.o : fastq-io.h fastq-io.c
	echo "Making fastq-io.o ..."
	$(CC) $(CFLAGS) fastq-io.c -c -lz -o fastq-io.o

fastq-dinuc-count : fastq-dinuc-count.c fastq-io.o
	echo "Making fastq-dinuc-count..."
	$(CC) $(CFLAGS) fastq-io.o fastq-dinuc-count.c -lz -o fastq-dinuc-count 
