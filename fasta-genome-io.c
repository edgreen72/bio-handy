#include "fasta-genome-io.h"

/* Takes filename as argument
   Returns true IFF filename ends in .gz
   False otherwise */
int is_gz( const char* fq_fn ) {
  size_t fn_len;
  fn_len = strlen( fq_fn );
  if ( (fq_fn[fn_len-3] == '.') &&
       (fq_fn[fn_len-2] == 'g') &&
       (fq_fn[fn_len-1] == 'z') ) {
    return 1;
  }
  return 0;
}

/* Takes filename
   Returns Fa_Src*; NULL => problem or nothing there
*/
Fa_Src* init_fasta_src( const char fn[] ) {
  Fa_Src* fa_source;
  if ( fn == NULL ) {
    return NULL;
  }
  fa_source = (Fa_Src*)malloc(sizeof( Fa_Src ));
  strcpy( fa_source->fn, fn );
  fa_source->n = 0;
  fa_source->seq_buffer = (char*)malloc(sizeof(char)*MAX_SEQ_LEN);
  fa_source->seq_buffer[0] = '\0';
  if ( is_gz( fn ) ) {
    fa_source->is_gz = 1;
    fa_source->fagz = gzopen( fa_source->fn, "r" );
    if ( fa_source->fagz == NULL ) {
      free( fa_source->seq_buffer );
      free( fa_source );
      return NULL;
    }
  }
  else {
    fa_source->is_gz = 0;
    fa_source->fafp = fileOpen( fa_source->fn, "r" );
    if ( fa_source->fafp == NULL ) {
      free( fa_source->seq_buffer );
      free( fa_source );
      return NULL;
    }
  }
  return fa_source;
}

/* get_next_fa
   Args: Fa_Src* fa_source - the source (uncompressed or gz) 
                 of some fasta data
   Returns: Seq* pointer to new sequence; NULL if there was
            any problem, like EOF
   Calls the right parser type (gz or regular) to read the next
   record and updates the fq_source->n if a fastq record is read
   correctly.
*/
Seq* get_next_fa( Fa_Src* fa_source, Genome* genome ) {
  Seq* seq;
  int status;
  seq = (Seq*)malloc(sizeof(Seq));
  if ( fa_source->is_gz ) {
    status = gzread_fasta( fa_source->fagz, seq,
			   fa_source->seq_buffer );

  }
  else {
    status = read_fasta( fa_source->fafp, seq,
			 fa_source->seq_buffer );
  }
  if ( status ) {
    free( seq );
    seq = NULL;
  }
  else {
    fa_source->n++;
    genome->seqs[ genome->n_seqs ] = seq;
    genome->n_seqs++;
  }
  return seq;
}

int close_fasta_src( Fa_Src* fa_source ) {
  if ( fa_source->is_gz ) {
    gzclose( fa_source->fagz );
  }
  else {
    fclose(fa_source->fafp);
  }
  return 0;
}

/* Args: FILE* fafp
         Seq* seq
         char* seq_buffer
   Returns: 0 - everything is copacetic
          non-zero if EOF or other problem
   Reads the next fastq sequence from a regular file handle
*/
int read_fasta( FILE* fafp, Seq* seq, char* seq_buffer ) {
  char c;
  size_t i = 0;
  c = fgetc(fafp);
  if ( c == '>' ) {
    c = fgetc(fafp);
    while ( !isspace(c) ) { // load up the ID
      seq->id[i++] = c;
      c = fgetc(fafp);
    }
    seq->id[i] = '\0';
    while( c != '\n' ) {
      c = fgetc(fafp);
    }
    i = 0;
    while( (c != '>') &&
           (c != EOF) &&
           (i < MAX_SEQ_LEN) ) {
      if ( isspace(c) ) {
        ;
      }
      else {
        c = toupper(c);
        seq_buffer[i++] = c;
      }
      c = fgetc( fafp );
    }
    seq_buffer[i] = '\0';
    ungetc(c, fafp);
  }
  else {
    if ( c == EOF ) {
      return -1;
    }
  }
  if ( i == MAX_SEQ_LEN ) {
    fprintf( stderr, "%s is truncated to %d\n", seq->id, MAX_SEQ_LEN );
  }
  /* Now, make space in seq for copying the sequence */
  seq->seq = (char*)malloc(sizeof(char)*(i+1));
  strcpy( seq->seq, seq_buffer );
  seq->len = i;
  return 0;
}

/* Args: gzFile gzfp
         Seq* seq
         char* seq_buffer
   Returns: 0 - everything is copacetic
          non-zero if EOF or other problem
   Reads the next fastq sequence from a gzFile handle
*/
int gzread_fasta( gzFile gzfp, Seq* seq, char* seq_buffer ) {
  char c;
  size_t i = 0;
  c = gzgetc(gzfp);
  if ( c == '>' ) {
    c = gzgetc(gzfp);
    while ( !isspace(c) ) { // load up the ID
      seq->id[i++] = c;
      c = gzgetc(gzfp);
    }
    seq->id[i] = '\0';
    while( c != '\n' ) {
      c = gzgetc(gzfp);
    }
    i = 0;
    while( (c != '>') &&
           (c != EOF) &&
           (i < MAX_SEQ_LEN) ) {
      if ( isspace(c) ) {
        ;
      }
      else {
        c = toupper(c);
        seq_buffer[i++] = c;
      }
      c = gzgetc( gzfp );
    }
    seq_buffer[i] = '\0';
    gzungetc(c, gzfp);
  }
  else {
    if ( c == EOF ) {
      return -1;
    }
  }
  if ( i == MAX_SEQ_LEN ) {
    fprintf( stderr, "%s is truncated to %d\n", seq->id, MAX_SEQ_LEN );
  }
  /* Now, make space in seq for copying the sequence */
  seq->seq = (char*)malloc(sizeof(char)*i);
  strcpy( seq->seq, seq_buffer );
  seq->len = i;
  return 0;
}

Seq* find_seq( Genome* genome, const char id[] ) {
  Seq** target;
  strcpy( genome->dummy->id, id );
  target = bsearch( &genome->dummy, genome->seqs, genome->n_seqs,
		    sizeof( Seq* ), chr_cmp );
  if ( target == NULL ) {
    return NULL;
  }
  else {
    return *target;
  }
}

int chr_cmp( const void *v1, const void *v2 ) {
  Seq** c1p = (Seq**) v1;
  Seq** c2p = (Seq**) v2;
  return strcmp( (*c1p)->id, (*c2p)->id );
}

Genome* init_genome( void ) {
  Genome* genome;
  genome = (Genome*)malloc(sizeof( Genome ));
  genome->seqs = (Seq**)malloc(sizeof(Seq*)*MAX_GENOME_SEQS);
  genome->dummy = (Seq*)malloc(sizeof(Seq));
  genome->n_seqs = 0;
  return genome;
}

/** fileOpen **/
FILE * fileOpen(const char *name, char access_mode[]) {
  FILE * f;
  f = fopen(name, access_mode);
  if (f == NULL) {
    fprintf( stderr, "%s\n", name);
    perror("Cannot open file");
    return NULL;
  }
  return f;
}

