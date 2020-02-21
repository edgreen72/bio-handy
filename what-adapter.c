#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>
#include <string.h>
#include <getopt.h>
#include "fastq-io.h"

#define DEBUG (0)
#define NUM_SEQ (100000)
#define ADAPT_LEN (65)
#define VERSION (3)

void help( char* adapter_root ) {
  printf( "what-adapter VERSION %d\n", VERSION );
  printf( "-f <fastq input file>\n" );
  printf( "-n <number of fastq sequences to examine; default = %d>\n", NUM_SEQ );
  printf( "-r <adapter root; default = %s>\n", adapter_root );
  printf( "-l <length of adapter to report; default = %d>\n", ADAPT_LEN );
  printf( "-v Verbose mode; report a bunch of stuff\n" );
  printf( "Report the adapter sequences seen in an input\n" );
  printf( "fastq file. The file may be gzipped or not.\n" );
  printf( "This program works by examining the first NUM_SEQ fastq\n" );
  printf( "records for the presence of the ADAPT_ROOT sequence.\n" );
  printf( "It then reports the most common full adapter sequence,\n" );
  printf( "i.e., the ADAPT_ROOT sequence and the most common sequence\n" );
  printf( "that follows it.\n" );
  printf( "This program is especially useful for situations where you\n" );
  printf( "have a short insert library and many reads have partial or\n" );
  printf( "adapter sequence.\n" );
  exit( 0 );
}

int main ( int argc, char* argv[] ) {
  extern char* optarg;
  char fq_in[MAX_FN_LEN+1] = {'\0'};
  int num_seq, adapt_len;
  char* adapter_root;
  char* adapter;
  FILE* fq;
  gzFile fqgz;
  FQ_Src* fq_source;
  FQ fq_seq;
  int ich;
  int verbose;
  
  /* Set up defaults */
  num_seq      = NUM_SEQ;
  adapt_len    = ADAPT_LEN;
  adapter_root = (char*)malloc(sizeof(char) * MAX_FQ_LEN);
  strcpy( adapter_root, "AGATCGGAAGAGC" );
  verbose = 0;

  /* Process input arguments */
  while( (ich=getopt( argc, argv, "f:n:l:r:v" )) != -1 ) {
    switch(ich) {
    case 'f' :
      strcpy( fq_in, optarg );
      break;
    case 'n' :
      num_seq = atoi( optarg );
      break;
    case 'l' :
      adapt_len = atoi( optarg );
      break;
    case 'r' :
      strcpy( adapter_root, optarg );
      break;
    case 'v' :
      verbose = 1;
      break;
    default :
      help(adapter_root);
    }
  }

  fq_source = init_fastq_src( fq_in );
  if ( fq_source == NULL ) {
    help(adapter_root);
  }
  while( get_next_fq( fq_source, &fq_seq ) == 0 ) {
    adapter = strstr( fq_seq.seq, adapter_root );
    if ( adapter != NULL ) {
      adapter[adapt_len] = '\0';
      printf( "%s\n", adapter );
    }
  }
  exit( 0 );
}

