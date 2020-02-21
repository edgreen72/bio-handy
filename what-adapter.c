#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>
#include <string.h>
#include <getopt.h>
#include "fastq-io.h"

#define DEBUG (0)
#define NUM_SEQ (100000)
#define L_DEF (167)

void help( void ) {
  printf( "what-adapter\n" );
  printf( "-f <fastq input file>\n" );
  printf( "-n <NUM_SEQ; default = %d>\n", NUM_SEQ );
  printf( "-r <ADAPT_ROOT; default = %s>\n", ADAPT_ROOT );
  printf( "-V verbose mode; tell me a bunch of other stuff\n" );

}

int main ( int argc, char* argv[] ) {
  extern char* optarg;
  char ADAPT_ROOT[] = "AGATCGGAAGAGC"
  char fq_in[MAX_FN_LEN+1] = {'\0'};
  //  char fq_fn[MAX_FN_LEN+1] = {'\0'};
  char* fq_fn;
  FILE* fq;
  gzFile fqgz;
  int length = L_DEF;
  DiNucArray DNA;
  FQ_Src* fq_source;
  FQ fq_seq;
  int ich;
  char delimiter = ':';
  
  while( (ich=getopt( argc, argv, "f:l:" )) != -1 ) {
    switch(ich) {
    case 'f' :
      strcpy( fq_in, optarg );
      break;
    case 'l' :
      length = atoi( optarg );
      break;
    default :
      help();
    }
  }

  DNA = init_DiNucArray( length );
  fq_fn = strtok( fq_in, &delimiter );
  fq_source = init_fastq_src( fq_fn );
  if ( fq_source == NULL ) {
    help();
  }
  while( fq_source != NULL ) {
    while( get_next_fq( fq_source, &fq_seq ) == 0 ) {
      if ( fq_seq.len == length ) {
	update_DNA( DNA, &fq_seq );
      }
    }
    fq_fn = strtok( NULL, &delimiter );
    fq_source = reset_fastq_src( fq_fn, fq_source );
  }
  write_DNA( DNA, length );
  exit( 0 );
}

void update_DNA( DiNucArray DNA, const FQ* fq_seq_p ) {
  size_t i = 0;
  size_t inx = 0;
  char b;
  
  for( i = 0; i < (fq_seq_p->len - 1); i++ ) {
    inx = get_dinuc_inx( &(fq_seq_p->seq[i]) );
    DNA->dnps[i]->dinuc_counts[inx]++;
  }
}

DiNucArray init_DiNucArray( const int length ) {
  DiNucArray DNA;
  size_t i, inx;

  DNA = (DiNucArray)malloc(sizeof(Dinuc_array));
  DNA->dnps = (DNP*)malloc(sizeof(DNP) * length);
  DNA->len = length;

  for( i = 0; i < length; i++ ) {
    DNA->dnps[i] = (DNP)malloc(sizeof(Dinucs));
    for( inx = 0; inx < 17; inx++ ) {
      DNA->dnps[i]->dinuc_counts[inx] = 0;
    }
  }
  return DNA;
}

void write_DNA( const DiNucArray DNA, const int length ) {
  size_t i, inx;
  unsigned int totals[17] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
			      0, 0, 0, 0, 0, 0, 0 };

  /* Find the total counts over all positions for each dinucleotide */
  for( i = 0; i < (DNA->len - 1); i++ ) {
    for( inx = 0; inx < 16; inx++ ) {
      totals[inx] += DNA->dnps[i]->dinuc_counts[inx];
    }
  }

  printf( "#Dinucleotides counts at each position on reads of length %d\n", length );
  printf( "#POS AA AC AG AT CA CC CG CT GA GC GG GT TA TC TG TT NN\n" );
  for( i = 0; i < (DNA->len - 1); i++ ) {
    printf( "%lu ", i );
    for( inx = 0; inx < 16; inx++ ) {
      printf( "%u ", DNA->dnps[i]->dinuc_counts[inx] );
    }
    printf( "%u\n", DNA->dnps[i]->dinuc_counts[inx] );
  }
  
  printf( "\n\n" );
  printf( "# Dinucleotides enrichment/depletion from average\n" );
  for( i = 0; i < (DNA->len - 1); i++ ) {
    printf( "%lu ", i );
    for( inx = 0; inx < 16; inx++ ) {
      if ( totals[inx] == 0 ) {
	printf( "0 " );
      }
      else {
	printf( "%.3f ", (float)((float)DNA->dnps[i]->dinuc_counts[inx] /
				 ((float)totals[inx] / ((float)length - 1.0))) );
      }
    }
    if ( totals[inx] == 0 ) {
      printf( "0\n" );
    }
    else {
      printf( "%.3f\n", (float)((float)DNA->dnps[i]->dinuc_counts[inx] /
				((float)totals[inx] / ((float)length - 1.0))) );
    }
  }
}


size_t get_dinuc_inx( const char* dinuc ) {
  char b1, b2;
  size_t inx = 0;
  b1 = toupper(dinuc[0]);
  b2 = toupper(dinuc[1]);

  switch( b1 ) {
  case 'A' :
    inx += 0;
    break;
  case 'C' :
    inx += 4;
    break;
  case 'G' :
    inx += 8;
    break;
  case 'T' :
    inx += 12;
    break;
  default :
    return 16;
  }

  switch(b2) {
  case 'A' :
    inx += 0;
    break;
  case 'C' :
    inx += 1;
    break;
  case 'G' :
    inx += 2;
    break;
  case 'T' :
    inx += 3;
    break;
  default :
    return 16;
  }
  return inx;  
}
