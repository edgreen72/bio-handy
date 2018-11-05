#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include <limits.h>
#include <zlib.h>
#define MAX_FN_LEN (2047)
#define MAX_ID_LEN (511)
#define MAX_SEQ_LEN (536870911)
#define MAX_GENOME_SEQS (1000000)

/* Data structures */
typedef struct seq {
  char id[MAX_ID_LEN + 1];
  char* seq;
  size_t len;
} Seq;

typedef struct genome {
  Seq** seqs;
  Seq* dummy; // used for searching
  size_t n_seqs;
} Genome;

typedef struct fa_src {
  char fn[MAX_FN_LEN+1];
  char* seq_buffer;
  int is_gz;
  gzFile fagz;
  FILE* fafp;
  size_t n;
} Fa_Src;

/* Function prototypes */
Genome* init_genome( void );
Fa_Src* init_fasta_src( const char fn[] );
Seq* get_next_fa( Fa_Src* fa_source, Genome* genome );
int read_fasta( FILE* fp, Seq* seq, char* seq_buffer );
int gzread_fasta( gzFile gzfp, Seq* seq, char* seq_buffer );
Seq* find_seq( Genome* genome, const char id[] );
int is_gz( const char* fn );
FILE* fileOpen( const char* name, char access_mode[] );
int close_fasta_src( Fa_Src* );
int chr_cmp( const void *v1, const void *v2 );
