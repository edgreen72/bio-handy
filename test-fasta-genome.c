#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>
#include <string.h>
#include <getopt.h>
#include "fasta-genome-io.h"

#define DEBUG (1)

void help( void ) {
  printf( "test-fasta-genome -f <fasta file> -I <ID of sequence to find>\n" );
  printf( "This program uses the fasta-genome-io code to parse an input\n" );
  printf( "fasta file representing a genome. It then attempts to find\n" );
  printf( "the sequence whose identifier is given via the -I option.\n" );
  exit( 0 );
}

int main( int argc, char* argv[] ) {
  extern char* optarg;
  int ich;
  char fa_in[ MAX_FN_LEN + 1 ]     = {'\0'};
  char target_id[ MAX_ID_LEN + 1 ] = {'\0'};
  Genome* genome;
  Seq* seq;
  Fa_Src* fa_src;

  while( (ich=getopt( argc, argv, "f:I:" ) ) != -1 ) {
    switch(ich) {
    case 'f' :
      strcpy( fa_in, optarg );
      break;
    case 'I' :
      strcpy( target_id, optarg );
      break;
    default :
      help();
    }
  }
  if ( (strlen( fa_in ) == 0) ||
       (strlen( target_id ) == 0) ) {
    help();
  }

  genome = init_genome();
  fa_src = init_fasta_src( fa_in );

  seq = get_next_fa( fa_src, genome );
  while( seq != NULL ) {
    if ( DEBUG ) {
      printf( "Saw %s length %lu\n", seq->id, seq->len );
    }
    seq = get_next_fa( fa_src, genome );
  }
  close_fasta_src( fa_src );
  qsort( genome->seqs, genome->n_seqs, sizeof(Seq*), chr_cmp );
  
  seq = find_seq( genome, target_id );
  if ( seq == NULL ) {
    printf( "Could not find %s in genome.\n",
	    target_id );
  }
  else {
    printf( "Found %s in genome. Length = %lu\n",
	    target_id, seq->len );
  }
  exit( 0 );
}
