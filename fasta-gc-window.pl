#!/usr/local/bin/perl

use Bio::SeqIO;
use Getopt::Std;
use strict;
use vars qw( $opt_f $opt_w );

my $bio_in = &init();
my ( $seq_o, $pos, $id, $len, $window_seq, $gc, $at, $window_seq );

while( $seq_o = $bio_in->next_seq() ) {
    $pos = 1;
    $id  = $seq_o->display_id();
    $len = $seq_o->length();
    while( ($pos + $opt_w) < $len ) {
	$window_seq = $seq_o->subseq( $pos, ($pos + $opt_w) );
	$gc = $window_seq =~ tr/GCgc//;
	$at = $window_seq =~ tr/ATat//;
	if ( ($gc + $at) > 0 ) {
	    printf( "%s %d %d\n", 
		    $id, 
		    $pos, 
		    int($gc/($gc+$at) * 100) );
	}
	$pos += $opt_w;
    }
    print( "\n\n" );
}

sub init {
    my $w_DEF = 100000;
    getopts( 'f:w:' );
    unless( -f $opt_f ) {
	printf( "fasta-gc-window.pl -f <fasta file> -w <window; DEF = $w_DEF>\n" );
	printf( "Makes a table of:\n" );
	printf( "1. Sequence ID\n" );
	printf( "2. Window start position\n" );
	printf( "3. %GC in window\n" );
	exit( 0 );
    }
    unless( defined( $opt_w ) ) {
	$opt_w = $w_DEF;
    }
    return Bio::SeqIO->new( '-file' => $opt_f,
			    '-format' => 'fasta' );
}
