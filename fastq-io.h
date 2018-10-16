#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include <limits.h>
#include <zlib.h>
#define MAX_FN_LEN (2047)
#define MAX_ID_LEN (511)
#define MAX_FQ_LEN (2047)

/* Data structures */
typedef struct fq {
  char id[ MAX_ID_LEN + 1];
  char seq[MAX_FQ_LEN + 1];
  char qual[MAX_FQ_LEN +1];
  size_t len;
} FQ;

typedef struct fqpair {
  FQ* fq1;
  FQ* fq2;
} FQPair;

typedef struct fq_src {
  char fn[MAX_FN_LEN + 1];
  int is_gz;
  gzFile fqgz;
  FILE* fqfp;
  size_t n; // number read so far
} FQ_Src;

typedef struct fqpair_src {
  FQ_Src* r1;
  FQ_Src* r2;
} FQPair_Src;

/* Function prototypes */
int is_gz( const char* fq_fn );
int get_next_fq( FQ_Src* fq_source, FQ* fq_seq );
int get_next_fqpair( FQPair_Src* fq_pair_source,
		     FQPair* fq_seq_pair );
FQ_Src* init_fastq_src( const char fn[] );
FQ_Src* reset_fastq_src( const char fn[], FQ_Src* fq_source );
int read_fastq( FILE* fp, FQ* fq_seq );
int gzread_fastq( gzFile gzfp, FQ* fq_seq );
FILE* fileOpen(const char* name, char access_mode[]);
