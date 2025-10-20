#include "kmer.h"

int seq2inx( char* seq, unsigned int k, unsigned int* inx ) {
  int i = 0;
  unsigned int inx_build = 0;
  while( i < k ) {
    inx_build = inx_build << 2;
    switch( seq[i] ) {
    case 'A' :
      inx_build += 0;
      break;
    case 'a' :
      inx_build += 0;
      break;
    case 'C' :
      inx_build += 1;
      break;
    case 'c' :
      inx_build += 1;
      break;
    case 'G' :
      inx_build += 2;
      break;
    case 'g' :
      inx_build += 2;
      break;
    case 'T' :
      inx_build += 3;
      break;
    case 't' :
      inx_build += 3;
      break;
    default :
      return 0;
    }
    i++;
  }
  *inx = inx_build;
  return 1;
}

int add_seq_to_KHA( KHA* kha, FQ* fq ) {
  unsigned int start_inx;
  unsigned int end_inx;
  unsigned int full_inx = 0;

  if ( seq2inx( fq->seq, kha->k, &start_inx ) == 0 ) {
    return 0;
  }
  if ( seq2inx( &fq->seq[fq->len - kha->k], kha->k, &end_inx )
       == 0 ) {
    return 0;
  }
  full_inx = start_inx;
  full_inx = full_inx << (2 * kha->k);
  full_inx += end_inx;
  
  /* First for this length? */
  if (kha->kaa[fq->len] == NULL ) {
    kha->kaa[fq->len] = init_kmer_array( kha->k );
  }
  kha->kaa[fq->len]->k_array[full_inx]++;
  return 1;
}

KA* init_kmer_array( const unsigned int k ) {
  size_t i;
  unsigned int array_size;
  KA* ka;
  ka = (KA*)malloc(sizeof( KA ));
  ka->k = k;
  array_size = 1 << (4*k);
  ka->k_array = (short unsigned int*)malloc(sizeof(short unsigned int)
					   * array_size);
  for( i = 0; i < array_size; i++ ) {
    ka->k_array[i] = 0;
  }
  return ka;
}

KHA* init_KHA( const unsigned int k ) {
  size_t i;
  KHA* kha;
  kha = (KHA*)malloc(sizeof(KHA));
  kha->k = k;
  kha->kaa_size = MAX_SEQ_LEN;
  kha->kaa = (KA**)malloc(sizeof(KA*) * (MAX_SEQ_LEN+1));
  for( i = 0; i <= MAX_SEQ_LEN; i++ ) {
    kha->kaa[i] = NULL;
  }
  return kha;
}
			  

/* 
 */
void output_kmer_table( KHA* kha ) {
  size_t len_i;
  size_t kmer_i;
  unsigned int i;
  unsigned int array_size;
  unsigned int count[MAX_SEQ_COUNT + 1];
  /* Zero out the count array that will be printed */
  for( i = 0; i <= MAX_SEQ_COUNT; i++ ) {
    count[i] = 0;
  }
  
  array_size = 1 << (4*kha->k);
  for( len_i = 0; len_i <= kha->kaa_size; len_i++ ) {
    /* Go through each kaa */
    if ( kha->kaa[len_i] != NULL ) {
      for( kmer_i = 0; kmer_i < array_size; kmer_i++ ) {
	if ( kha->kaa[len_i]->k_array[kmer_i] > MAX_SEQ_COUNT ) {
	  count[MAX_SEQ_COUNT]++;
	}
	else {
	  count[kha->kaa[len_i]->k_array[kmer_i]]++;
	}
      }
    }
  }
  for( i = 1; i <= MAX_SEQ_COUNT; i++ ) {
    printf( "%u %u\n", i, count[i] );
  }
}
