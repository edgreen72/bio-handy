#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>
#include <string.h>
#include <getopt.h>
#include "fastq-io.h"

#define DEBUG (0)
#define L_DEF (167)

typedef struct dinucs {
  unsigned int dinuc_counts[17];
} Dinucs;
typedef struct dinucs* DNP;

typedef struct dinuc_array {
  DNP* dnps;
  size_t len;
} Dinuc_array;
typedef struct dinuc_array* DiNucArray;

DiNucArray init_DiNucArray( const int length );
void update_DNA( DiNucArray DNA, const FQ* fq_seq_p );
size_t get_dinuc_inx( const char* dinuc );
void write_DNA( const DiNucArray DNA, const int length );


void help( void ) {
  printf( "fastq-dinuc-count -f <fastq file> -l <length>\n" );
  printf( "Makes a table of the observed dinucleotides in sequences\n" );
  printf( "of a defined length in a fastq file. Input file can be\n" );
  printf( "gzipped or not.\n" );
  printf( "Output is of this format:\n" );
  printf( "Position AA AC AG AT CA CC CG CT GA GC GG GT TA TC TG TT NN\n" );
  printf( "Position is the 0-indexed starting position of the dinucleotide.\n" );
  printf( "The NN column contains the count of all dinucleotides at the\n" );
  printf( "given position that cannot be counted in another category, i.e.,\n" );
  printf( "any dinucleotide composed of anything other that two\n" );
  printf( "A, C, G, or T bases.\n" );
  exit( 0 );
}

int main ( int argc, char* argv[] ) {
  extern char* optarg;
  char fq_fn[MAX_FN_LEN+1] = {'\0'};
  FILE* fq;
  gzFile fqgz;
  int length = L_DEF;
  DiNucArray DNA;
  FQ_Src* fq_source;
  FQ fq_seq;
  int ich;
  
  while( (ich=getopt( argc, argv, "f:l:" )) != -1 ) {
    switch(ich) {
    case 'f' :
      strcpy( fq_fn, optarg );
      break;
    case 'l' :
      length = atoi( optarg );
      break;
    default :
      help();
    }
  }

  DNA = init_DiNucArray( length );
  fq_source = init_fastq_src( fq_fn );

  if ( fq_source == NULL ) {
    help();
  }
  
  while( get_next_fq( fq_source, &fq_seq ) == 0 ) {
    if ( fq_seq.len == length ) {
      update_DNA( DNA, &fq_seq );
    }
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
  printf( "#Dinucleotides counts at each position on reads of length %d\n", length );
  printf( "#POS AA AC AG AT CA CC CG CT GA GC GG GT TA TC TG TT NN\n" );
  for( i = 0; i < (DNA->len - 1); i++ ) {
    printf( "%lu ", i );
    for( inx = 0; inx < 16; inx++ ) {
      printf( "%u ", DNA->dnps[i]->dinuc_counts[inx] );
    }
    printf( "%u\n", DNA->dnps[i]->dinuc_counts[inx] );
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
