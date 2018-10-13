#include "fastq-io.h"

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


/* get_next_fqpair
   Args: FQPair_Src* fq_pair_source - the input sources for a fastq read pair
         FQPair* fq_seq_pair - the place to put the read pair data
   Returns: 0 => data read; everything copacetic
           -1 => EOF or other problem; stop trying
 */
int get_next_fqpair( FQPair_Src* fq_pair_source, FQPair* fq_seq_pair ) {
  if ( (get_next_fq( fq_pair_source->r1, fq_seq_pair->fq1 ) == 0) &&
       (get_next_fq( fq_pair_source->r2, fq_seq_pair->fq2 ) == 0) ) {
    return 0;
  }
  return -1;
}

FQ_Src* init_fastq_src( const char fn[] ) {
  FQ_Src* fq_source;
  fq_source = (FQ_Src*)malloc(sizeof( FQ_Src ));
  strcpy( fq_source->fn, fn );
  fq_source->n = 0;
  if ( is_gz( fn ) ) {
    fq_source->is_gz = 1;
    fq_source->fqgz = gzopen( fq_source->fn, "r" );
  }
  else {
    fq_source->is_gz = 0;
    fq_source->fqfp = fileOpen( fq_source->fn, "r" );
  }
  return fq_source;
}

/* get_next_fq
   Args: FQ_Src* fq_source - the source (uncompressed or gz) of some
                             fastq data
         FQ* fq_seq - the place to put the data
   Returns: 0 => read next fastq record; everything copacetic
           -1 => EOF or other problem; stop trying on this source
   Calls the right parser type (gz or regular) to read the next
   record and updates the fq_source->n if a fastq record is read
   correctly.
*/
int get_next_fq( FQ_Src* fq_source, FQ* fq_seq ) {
  if ( fq_source->is_gz ) {
    if ( gzread_fastq( fq_source->fqgz, fq_seq ) ) {
      return -1;
    }
    else {
      fq_source->n++;
      return 0;
    }
  }
  else {
    if ( read_fastq( fq_source->fqfp, fq_seq ) ) {
      return -1;
    }
    else {
      fq_source->n++;
      return 0;
    }
  }
}

/* Args: FQ_Src* fq_source - contains the filepointer to fastq file
         FQ* fq_seq - data structure to put the fastq data for the next fastq seq
   Returns: 0 - everything is copacetic
          non-zero if EOF or other problem
   Reads the next fastq sequence from a regular file handle
*/
int read_fastq( FILE* fastq, FQ* fq_seq ) {
  char c;
  size_t i;
  
  c = fgetc( fastq );
  if ( c == EOF ) return -1;
  if ( c != '@' ) {
    fprintf( stderr, "fastq record not beginning with @\n" );
    return -1;
  }

  /* get identifier */
  i = 0;
  while( (!isspace(c=fgetc( fastq ) ) &&
          (i < MAX_ID_LEN) ) ) {
    if ( c == EOF ) {
      return -1;
    }
    fq_seq->id[i] = c;
    i++;
    if ( i == MAX_ID_LEN ) {
      /* Id is too long - truncate it now */
      fq_seq->id[i] = '\0';
    }
  }
  fq_seq->id[i] = '\0';

  /* Now, everything else on the line is description (if anything)
     although fastq does not appear to formally support description */
  while ( (c != '\n') &&
          (c != EOF) ) {
    c = fgetc( fastq );
  }
  i = 0;

  /* Now, read the sequence. This should all be on a single line */
  i = 0;
  c = fgetc( fastq );
  while ( (c != '\n') &&
          (c != EOF) &&
          (i < MAX_FQ_LEN) ) {
    if ( isspace(c) ) {
      ;
    }
    else {
      c = toupper( c );
      fq_seq->seq[i++] = c;
    }
    c = fgetc( fastq );
  }
  fq_seq->seq[i] = '\0';
  fq_seq->len = i;
  
  /* If the reading stopped because the sequence was longer than
     MAX_FQ_LEN, then we need to advance the file pointer
     past this line */
  if ( i == MAX_FQ_LEN ) {
    while ( (c != '\n') &&
            (c != EOF) ) {
      c = fgetc( fastq );
    }
  }

  /* Now, read the quality score header */
  c = fgetc( fastq );
  if ( c != '+' ) {
    fprintf( stderr, "Problem reading quality line for %s\n", fq_seq->id );
    return 1;
  }
  /* Zip through the rest of the line, it should be the same identifier
     as before or blank */
  c = fgetc( fastq );
  while( (c != '\n') &&
         (c != EOF) ) {
    c = fgetc( fastq );
  }

  /* Now, get the quality score line */
  c = fgetc( fastq );
  i = 0;
  while( (c != '\n') &&
         (c != EOF) &&
         (i < MAX_FQ_LEN) ) {
    if ( isspace(c) ) {
      ;
    }
    else {
      fq_seq->qual[i++] = c;
    }
    c = fgetc( fastq );
  }
  fq_seq->qual[i] = '\0';

  /* If the reading stopped because the sequence was longer than
     MAX_FQ_LEN, then we need to advance the file pointer
     past this line */
  if ( i == MAX_FQ_LEN ) {
    while ( (c != '\n') &&
            (c != EOF) ) {
      c = fgetc( fastq );
    }
  }

  if ( c == EOF ) {
    return -1;
  }
  return 0;
}

/* Args: gzFile fastq - file pointer to a fastq file ready to read 
                       the next sequence
         FQ* fq_seq - data structure to put the fastq data for the next fastq seq
   Returns: 0 - everything is copacetic
          non-zero if EOF or other problem
   Reads the next fastq sequence from a gzFile (gzipped) filehandle
*/
int gzread_fastq( gzFile fastq, FQ* fq_seq ) {
		  //char id[], char seq[], char qual[], size_t* len ) {
  char c;
  size_t i;
  c = gzgetc( fastq );
  if ( c == EOF ) return -1;
  if ( c != '@' ) {
    fprintf( stderr, "fastq record not beginning with @\n" );
    return -1;
  }

  /* get identifier */
  i = 0;
  while( (!isspace(c=gzgetc( fastq ) ) &&
          (i < MAX_ID_LEN) ) ) {
    if ( c == EOF ) {
      return -1;
    }
    fq_seq->id[i] = c;
    i++;
    if ( i == MAX_ID_LEN ) {
      /* Id is too long - truncate it now */
      fq_seq->id[i] = '\0';
    }
  }
  fq_seq->id[i] = '\0';

  /* Now, everything else on the line is description (if anything)
     although fastq does not appear to formally support description */
  while ( (c != '\n') &&
          (c != EOF) ) {
    c = gzgetc( fastq );
  }
  i = 0;

  /* Now, read the sequence. This should all be on a single line */
  i = 0;
  c = gzgetc( fastq );
  while ( (c != '\n') &&
          (c != EOF) &&
          (i < MAX_FQ_LEN) ) {
    if ( isspace(c) ) {
      ;
    }
    else {
      c = toupper( c );
      fq_seq->seq[i++] = c;
    }
    c = gzgetc( fastq );
  }
  fq_seq->seq[i] = '\0';
  fq_seq->len = i;
  
  /* If the reading stopped because the sequence was longer than
     MAX_FQ_LEN, then we need to advance the file pointer
     past this line */
  if ( i == MAX_FQ_LEN ) {
    while ( (c != '\n') &&
            (c != EOF) ) {
      c = gzgetc( fastq );
    }
  }

  /* Now, read the quality score header */
  c = gzgetc( fastq );
  if ( c != '+' ) {
    fprintf( stderr, "Problem reading quality line for %s\n", fq_seq->id );
    return 1;
  }
  /* Zip through the rest of the line, it should be the same identifier
     as before or blank */
  c = gzgetc( fastq );
  while( (c != '\n') &&
         (c != EOF) ) {
    c = gzgetc( fastq );
  }

  /* Now, get the quality score line */
  c = gzgetc( fastq );
  i = 0;
  while( (c != '\n') &&
         (c != EOF) &&
         (i < MAX_FQ_LEN) ) {
    if ( isspace(c) ) {
      ;
    }
    else {
      fq_seq->qual[i++] = c;
    }
    c = gzgetc( fastq );
  }
  fq_seq->qual[i] = '\0';

  /* If the reading stopped because the sequence was longer than
     MAX_FQ_LEN, then we need to advance the file pointer
     past this line */
  if ( i == MAX_FQ_LEN ) {
    while ( (c != '\n') &&
            (c != EOF) ) {
      c = gzgetc( fastq );
    }
  }

  if ( c == EOF ) {
    return -1;
  }
  return 0;
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
