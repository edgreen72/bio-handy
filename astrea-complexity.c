#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>
#include <string.h>
#include <getopt.h>
#include "fastq-io.h"
#include "kmer.h"
#define VERSION (1)
#define DEF_KMER_LEN (6)
#define MAX_TO_READ (1000000)

void help( void );

int main ( int argc, char* argv[] ) {
  extern char* optarg;
  char fq_fn[MAX_FN_LEN + 1];
  KHA* kha;
  FILE* fq;
  gzFile fqgz;
  FQ_Src* fq_source;
  FQ fq_seq;
  int i;
  int ich;
  unsigned int k           = DEF_KMER_LEN;
  unsigned int total       = 0;
  unsigned int unreadable  = 0;
  unsigned int max_to_read = MAX_TO_READ;
  if ( argc == 1 ) {
    help();
  }

  while( (ich=getopt( argc, argv, "f:k:m:" )) != -1 ) {
    switch(ich) {
    case 'f' :
      strcpy( fq_fn, optarg );
      break;
    case 'k' :
      k = atoi( optarg );
      break;
    case 'm' :
      max_to_read = atoi( optarg );
      break;
    default :
      help();
    }
  }

  fq_source = init_fastq_src( fq_fn );
  if ( fq_source == NULL ) {
    help();
  }

  kha = init_KHA( k );
  if ( kha == NULL ) {
    help();
  }
  
  while( (get_next_fq( fq_source, &fq_seq ) == 0) &&
	 (total < max_to_read) ) {
    if ( add_seq_to_KHA( kha, &fq_seq ) ) {
      total++;
    }
    else {
      unreadable++;
    }
  }


  printf( "# Complexity analysis of %s\n", fq_fn );
  printf( "# Total sequences with start and end kmers: %u\n",
	  total );
  printf( "# Total sequences without start or end kmers: %u\n",
	  unreadable );
  output_kmer_table( kha );
  exit( 0 );
}

void help( void ) {
  printf("astrea-complexity V %d\n", VERSION );
  printf("    -f <fastq file>\n" );
  printf("    -k <kmer length; default = %d", DEF_KMER_LEN);
  printf("    -m <max sequences to examine; def = %d\n",
	 MAX_TO_READ );
  printf("Makes a histogram of how many sequences are seen\n" );
  printf("each specific number of times.\n" );
  printf("A sequence is the same if its length is the same\n" );
  printf("and the first k and last k bases of the sequence\n" );
  printf("are the same.\n");
  printf("This is meant to be run on merged sequence data,\n" );
  printf("i.e., not read pairs.\n" );
  exit( 0 );
  
} 

