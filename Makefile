CC=gcc
CFLAGS=-gdwarf-2 -g

fastq-io.o : fastq-io.h fastq-io.c
	echo "Making fastq-io.o ..."
	$(CC) $(CFLAGS) fastq-io.c -c -lz -o fastq-io.o
