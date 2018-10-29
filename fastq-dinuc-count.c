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
void write_DNA( const DiNucArray DNA, const int length, const char out_fn_root[] );
void make_gnuplot_plot( const char out_fn_root[], const int length );

void help( void ) {
  printf( "fastq-dinuc-count -f <fastq file> -l <length>\n" );
  printf( "                  -e <write output files to this name>\n" );
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
  printf( "NOTES: -f option can be a colon delimited list of fastq filenames\n" );
  printf( "       If -e option is given, then a .dat and .eps filename are\n" );
  printf( "       made with this prefix. These files contain the data table (.dat)\n" );
  printf( "       and an Encapsulated Postscript File (.eps) with an image of\n" );
  printf( "       the results.\n" );
  exit( 0 );
}

int main ( int argc, char* argv[] ) {
  extern char* optarg;
  char fq_in[MAX_FN_LEN+1] = {'\0'};
  char out_fn_root[MAX_FN_LEN+1] = {'\0'};
  char* fq_fn;
  FILE* fq;
  gzFile fqgz;
  int length = L_DEF;
  DiNucArray DNA;
  FQ_Src* fq_source;
  FQ fq_seq;
  int ich;
  int make_plot = 0;
  char delimiter = ':';
  
  while( (ich=getopt( argc, argv, "f:l:e:" )) != -1 ) {
    switch(ich) {
    case 'f' :
      strcpy( fq_in, optarg );
      break;
    case 'l' :
      length = atoi( optarg );
      break;
    case 'e' :
      strcpy( out_fn_root, optarg );
      make_plot = 1;
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
    fprintf( stderr, "Examining %s...\n", fq_source->fn );
    while( get_next_fq( fq_source, &fq_seq ) == 0 ) {
      if ( fq_seq.len == length ) {
	update_DNA( DNA, &fq_seq );
      }
    }
    fq_fn = strtok( NULL, &delimiter );
    fq_source = reset_fastq_src( fq_fn, fq_source );
  }
  write_DNA( DNA, length, out_fn_root );
  if ( make_plot ) {
    make_gnuplot_plot( out_fn_root, length );
  }
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

void write_DNA( const DiNucArray DNA, const int length, const char out_fn_root[] ) {
  size_t i, inx;
  FILE* out_fh;
  char dat_fn[ MAX_FN_LEN+1 ] = {'\0'};
  unsigned int totals[17] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
			      0, 0, 0, 0, 0, 0, 0 };

  if ( strlen( out_fn_root ) == 0 ) {
    out_fh = stdout;
  }
  else {
    strcpy( dat_fn, out_fn_root );
    strcat( dat_fn, ".dat" );
    out_fh = fileOpen( dat_fn, "w" );
    if( out_fh == NULL ) {
      fprintf( stderr, "Cannot write to %s\n", dat_fn );
      exit( 1 );
    }
  }
  
  /* Find the total counts over all positions for each dinucleotide */
  for( i = 0; i < (DNA->len - 1); i++ ) {
    for( inx = 0; inx < 16; inx++ ) {
      totals[inx] += DNA->dnps[i]->dinuc_counts[inx];
    }
  }

  fprintf( out_fh, "#Dinucleotides counts at each position on reads of length %d\n", length );
  fprintf( out_fh, "#POS AA AC AG AT CA CC CG CT GA GC GG GT TA TC TG TT NN\n" );
  for( i = 0; i < (DNA->len - 1); i++ ) {
    fprintf( out_fh, "%lu ", i );
    for( inx = 0; inx < 16; inx++ ) {
      fprintf( out_fh, "%u ", DNA->dnps[i]->dinuc_counts[inx] );
    }
    fprintf( out_fh, "%u\n", DNA->dnps[i]->dinuc_counts[inx] );
  }
  
  fprintf( out_fh, "\n\n" );
  fprintf( out_fh, "# Dinucleotides enrichment/depletion from average\n" );
  for( i = 0; i < (DNA->len - 1); i++ ) {
    fprintf( out_fh, "%lu ", i );
    for( inx = 0; inx < 16; inx++ ) {
      if ( totals[inx] == 0 ) {
	fprintf( out_fh, "0 " );
      }
      else {
	fprintf( out_fh, "%.3f ", (float)((float)DNA->dnps[i]->dinuc_counts[inx] /
					  ((float)totals[inx] / ((float)length - 1.0))) );
      }
    }
    if ( totals[inx] == 0 ) {
      fprintf( out_fh, "0\n" );
    }
    else {
      fprintf( out_fh, "%.3f\n", (float)((float)DNA->dnps[i]->dinuc_counts[inx] /
					 ((float)totals[inx] / ((float)length - 1.0))) );
    }
  }
  fclose( out_fh );
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

void make_gnuplot_plot( const char out_fn_root[], const int length ) {
  char eps_fn[MAX_FN_LEN+1] = {'\0'};
  char dat_fn[MAX_FN_LEN+1] = {'\0'};
  size_t num_commands = 28;
  size_t i;
  FILE* gnuplotPipe;
  char* commands[] = {
    "set terminal postscript eps color solid \"Arial\" 14;",
    "set style line 1 lt 1 lw 2 pt 7 lc rgb \"#483D8B\"; #darkslateblue",
    "set style line 2 lt 1 lw 2 pt 7 lc rgb \"#FF8C00\"; #aqua",
    "set style line 3 lt 1 lw 2 pt 7 lc rgb \"#FFD700\"; # gold",
    "set style line 4 lt 1 lw 2 pt 7 lc rgb \"#1E90FF\"; # dodgerblue",
    "set grid;",
    "set key inside center top horizontal;",
    "set xtic 1;",
    "set ytic 1;",
    "set multiplot;",
    "set size 0.5, 0.25;",
    "set tmargin 0.95;",
    "set origin 0, 0.75;",
    "plot [0:15][-1.6:1.6] datafile index 1 using ($1+1):(log($2)) t \"AA\" w lp ls 1, '' index 1 using ($1+1):(log($3)) t \"AC\" w lp ls 2, '' index 1 using ($1+1):(log($4)) t \"AG\" w lp ls 3, '' index 1 using ($1+1):(log($5)) t \"AT\" w lp ls 4;",

    "set origin 0.5, 0.75;",
    "plot [length-15:length][-1.6:1.6] datafile index 1 using ($1+1):(log($2)) t \"AA\" w lp ls 1, '' index 1 using ($1+1):(log($3)) t \"AC\" w lp ls 2, '' index 1 using ($1+1):(log($4)) t \"AG\" w lp ls 3, '' index 1 using ($1+1):(log($5)) t \"AT\" w lp ls 4;",

    "set origin 0, 0.5;",
    "plot [0:15][-1.6:1.6] datafile index 1 using ($1+1):(log($6)) t \"CA\" w lp ls 1, '' index 1 using ($1+1):(log($7)) t \"CC\" w lp ls 2, '' index 1 using ($1+1):(log($8)) t \"CG\" w lp ls 3, '' index 1 using ($1+1):(log($9)) t \"CT\" w lp ls 4;",

    "set origin 0.5, 0.5;",
    "plot [length-15:length][-1.6:1.6] datafile index 1 using ($1+1):(log($6)) t \"CA\" w lp ls 1, '' index 1 using ($1+1):(log($7)) t \"CC\" w lp ls 2, '' index 1 using ($1+1):(log($8)) t \"CG\" w lp ls 3, '' index 1 using ($1+1):(log($9)) t \"CT\" w lp ls 4;",

    "set origin 0, 0.25;",
    "plot [0:15][-1.6:1.6] datafile index 1 using ($1+1):(log($10)) t \"GA\" w lp ls 1, '' index 1 using ($1+1):(log($11)) t \"GC\" w lp ls 2, '' index 1 using ($1+1):(log($12)) t \"GG\" w lp ls 3, '' index 1 using ($1+1):(log($14)) t \"GT\" w lp ls 4;",

    "set origin 0.5, 0.25;",
    "plot [length-15:length][-1.6:1.6] datafile index 1 using ($1+1):(log($10)) t \"GA\" w lp ls 1, '' index 1 using ($1+1):(log($11)) t \"GC\" w lp ls 2, '' index 1 using ($1+1):(log($12)) t \"GG\" w lp ls 3, '' index 1 using ($1+1):(log($13)) t \"GT\" w lp ls 4;",

    "set origin 0, 0;",
    "plot [0:15][-1.6:1.6] datafile index 1 using ($1+1):(log($14)) t \"TA\" w lp ls 1, '' index 1 using ($1+1):(log($15)) t \"TC\" w lp ls 2, '' index 1 using ($1+1):(log($16)) t \"TG\" w lp ls 3, '' index 1 using ($1+1):(log($17)) t \"TT\" w lp ls 4;",

    "set origin 0.5, 0;",
    "plot [length-15:length][-1.6:1.6] datafile index 1 using ($1+1):(log($14)) t \"TA\" w lp ls 1, '' index 1 using ($1+1):(log($15)) t \"TC\" w lp ls 2, '' index 1 using ($1+1):(log($16)) t \"TG\" w lp ls 3, '' index 1 using ($1+1):(log($17)) t \"TT\" w lp ls 4;",

  };

  /* Set up filename strings for input & output files */
  strcpy( eps_fn, out_fn_root );
  strcpy( dat_fn, out_fn_root );
  strcat( eps_fn, ".eps" );
  strcat( dat_fn, ".dat" );

  /* Open filehandle to give gnuplot commands */
  gnuplotPipe = popen( "gnuplot", "w" );
  
  /* Print some commands to gnuplot */
  fprintf( gnuplotPipe, "set output \"%s\";\n", eps_fn );
  fprintf( gnuplotPipe, "datafile = \"%s\";\n", dat_fn );
  fprintf( gnuplotPipe, "length = %d;\n", length );

  for( i = 0; i < num_commands; i++ ) {
    fprintf( gnuplotPipe, "%s \n", commands[i] );
  }
  
  /* Flush & close */
  fflush( gnuplotPipe );
  pclose( gnuplotPipe );
}

/*
set terminal postscript eps color solid "Arial" 14;
set output "t5.eps";
datafile = "t5.dat";
set style line 1 lt 1 lw 2 pt 7 lc rgb "#483D8B"; #darkslateblue
set style line 2 lt 1 lw 2 pt 7 lc rgb "#FF8C00"; #aqua
set style line 3 lt 1 lw 2 pt 7 lc rgb "#FFD700"; # gold
set style line 4 lt 1 lw 2 pt 7 lc rgb "#1E90FF"; # dodgerblue

set grid;
set key inside center top horizontal;
set xtic 1;
set ytic 1;
set multiplot;
set size 0.5, 0.25;
set tmargin 0.95;
#set bmargin 0.99;
set origin 0, 0.75;
plot [0:15][-1.6:1.6] datafile index 1 using ($1+1):(log($2)) t "AA" w lp ls 1,\
     '' index 1 using ($1+1):(log($3)) t "AC" w lp ls 2,\
     '' index 1 using ($1+1):(log($4)) t "AG" w lp ls 3,\
     '' index 1 using ($1+1):(log($5)) t "AT" w lp ls 4;

set origin 0.5, 0.75;
plot [26:41][-1.6:1.6] datafile index 1 using ($1+1):(log($2)) t "AA" w lp ls 1,\
     '' index 1 using ($1+1):(log($3)) t "AC" w lp ls 2,\
     '' index 1 using ($1+1):(log($4)) t "AG" w lp ls 3,\
     '' index 1 using ($1+1):(log($5)) t "AT" w lp ls 4;

set origin 0, 0.5;
plot [0:15][-1.6:1.6] datafile index 1 using ($1+1):(log($6)) t "CA" w lp ls 1,\
     '' index 1 using ($1+1):(log($7)) t "CC" w lp ls 2,\
     '' index 1 using ($1+1):(log($8)) t "CG" w lp ls 3,\
     '' index 1 using ($1+1):(log($9)) t "CT" w lp ls 4;

set origin 0.5, 0.5;
plot [26:41][-1.6:1.6] datafile index 1 using ($1+1):(log($6)) t "CA" w lp ls 1,\
     '' index 1 using ($1+1):(log($7)) t "CC" w lp ls 2,\
     '' index 1 using ($1+1):(log($8)) t "CG" w lp ls 3,\
     '' index 1 using ($1+1):(log($9)) t "CT" w lp ls 4;

set origin 0, 0.25;
plot [0:15][-1.6:1.6] datafile index 1 using ($1+1):(log($10)) t "GA" w lp ls 1,\
     '' index 1 using ($1+1):(log($11)) t "GC" w lp ls 2,\
     '' index 1 using ($1+1):(log($12)) t "GG" w lp ls 3,\
     '' index 1 using ($1+1):(log($13)) t "GT" w lp ls 4;

set origin 0.5, 0.25;
plot [26:41][-1.6:1.6] datafile index 1 using ($1+1):(log($10)) t "GA" w lp ls 1,\
     '' index 1 using ($1+1):(log($11)) t "GC" w lp ls 2,\
     '' index 1 using ($1+1):(log($12)) t "GG" w lp ls 3,\
     '' index 1 using ($1+1):(log($13)) t "GT" w lp ls 4;

set origin 0, 0;
plot [0:15][-1.6:1.6] datafile index 1 using ($1+1):(log($14)) t "TA" w lp ls 1,\
     '' index 1 using ($1+1):(log($15)) t "TC" w lp ls 2,\
     '' index 1 using ($1+1):(log($16)) t "TG" w lp ls 3,\
     '' index 1 using ($1+1):(log($17)) t "TT" w lp ls 4;

set origin 0.5, 0;
plot [26:41][-1.6:1.6] datafile index 1 using ($1+1):(log($14)) t "TA" w lp ls 1,\
     '' index 1 using ($1+1):(log($15)) t "TC" w lp ls 2,\
     '' index 1 using ($1+1):(log($16)) t "TG" w lp ls 3,\
     '' index 1 using ($1+1):(log($17)) t "TT" w lp ls 4;
*/
