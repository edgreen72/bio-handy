#ifndef KMER
#define KMER

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "fastq-io.h"
#define MAX_SEQ_COUNT (256)
#define MAX_SEQ_LEN (511)

/* KA is the array that will keep counts of each k-mer of length k.
   The array is indexed by converting the kmer to a number using
   the formula A=00, C=01, G=10, T=11. The kmer is converted to
   a bitstring. That bitstring is interpreted as an unsigned int.
 */
typedef struct kmer_array {
  unsigned int k;
  short unsigned int* k_array;
} KA;

/* KHA is an array with pointers so KA
   This array can be indexed by the length of the sequence or
   any arbitrary thing.
 */
typedef struct kmer_array_array {
  unsigned int k;
  unsigned int kaa_size; // length of KA** kaa
  unsigned int n_entries; // how many seqs in KA** kaa
  KA** kaa;
} KHA;

/* Takes the sequence, k-mer length, and the index to compute
   Converts the first k bases of the seq to bitstring/unsigned int
   using the formula A=00, C=01, G=10, T=11.
   returns 1 if copacetic.
   Returns 0 if there is a non ACGT character or sequence is
   less than k in length
 */
int seq2inx( char* seq, unsigned int k, unsigned int* inx );

/* Takes the KHA* and FQ* seq
   Increments the correct length and k-mer for this sequence
   given the k value of KHA*. Uses the beginning and ending
   kmer */
int add_seq_to_KHA( KHA* kha, FQ* fq );

/* Initializes the KHA, returns pointer to it */
KHA* init_KHA( const unsigned int k );
/* init_kmer_array 
   Returns pointer for kmer array of size 2*k
*/
KA* init_kmer_array( const unsigned int k );

/* 
 */
void output_kmer_table( KHA* kha );

#endif
