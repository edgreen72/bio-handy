#!/usr/bin/perl

use Bio::DB::Sam;
use Getopt::Std;
use vars qw( $opt_f $opt_b $opt_r $opt_F $opt_R );
use strict;

my ( $sam, $feature, $ref, $matches, $query, $i, $rcref, $rcquery );
my ( @ref, @query, @rref, @rquery );
my (@SUBS, @RSUBS );
my @BASES = ('A', 'C', 'G', 'T' );
my %RC = ('A' => 'T', 'C' => 'G', 'G' => 'C', 'T' => 'A',
	  'a' => 't', 'c' => 'g', 'g' => 'c', 't' => 'a',
	  'N' => 'N', 'n' => 'n', 'X' => 'x', '-' => '-');
	  
&init();

my $sam = Bio::DB::Sam->new( -bam   => $opt_b,
			     -fasta => $opt_f );

foreach $feature ( $sam->features('match') ) {
    unless( $feature->unmapped ) {
	($ref, $matches, $query) = $feature->padded_alignment();
	if ( $feature->strand == -1 ) {
	    $ref = &revcom( $ref );
	    $query = &revcom( $query );
	}
	@ref   = split( '', $ref );
	@query = split( '', $query );
	for( $i = 0; $i < $opt_r; $i++ ) {
	    $SUBS[$i]->{$ref[$i]}->{$query[$i]}++;
	}

	### Now the other way
#	$rcref   = &revcom( $ref );
#	$rcquery = &revcom( $query );
	####    NOTE: NOT reverse complement. We want to see what's on the top strand
	####          here, too - not the reverse complement of it, which would be
	####          the bottom strand
	@rref = reverse (split( '', $ref ));
	@rquery = reverse (split( '', $query ));
	for( $i = 0; $i < $opt_r; $i++ ) {
	    $RSUBS[$i]->{$rref[$i]}->{$rquery[$i]}++;
	}
    }
}

&output( '### Forward read substitution counts',
	 \@SUBS );

&output( "\n\n### Reverse read substitution counts",
	 \@RSUBS );

0;

sub output {
    my $header_string = shift;
    my $s_p = shift;
    my ( $i, $ref_base, $query_base, $count );
    print( "$header_string\n" );

    for( $i = 0; $i < $opt_r; $i++ ) {
	print( "$i" );
	foreach $ref_base (@BASES) {
	    foreach $query_base (@BASES) {
		if ( $s_p->[$i]->{$ref_base}->{$query_base} ) {
		    $count = $s_p->[$i]->{$ref_base}->{$query_base};
		}
		else {
		    $count = 0;
		}
		printf(" $count");
	    }
	}
	print( "\n" );
    }
}

#foreach $target ($sam->seq_ids) {
#    $segment = $sam->segment( -seq_id => $target );
#    my $iterator = $segment->features( -iterator => 1 );
#    while( my $align = $iterator->next_seq ) {
#	print "yes";
#    }
#}

sub revcom {
    my $seq = shift;
    my @seq = split( '', $seq );
    my @rcseq = map { $RC{$_} } reverse @seq;

    return join( '', @rcseq );
}

sub init {
    my $r_DEF = 15;
    getopts( 'f:b:r:' );
    unless( -f $opt_f &&
	    -f $opt_b ) {
	print( "pss-bam.pl -f <fasta file> -b <bam file>\n" );
	print( "           -r <region length; default = $r_DEF\n" );
	print( "Tests the Bio::DB::Sam module\n" );
	exit( 0 );
    }
    unless( defined( $opt_r ) ) {
	$opt_r = $r_DEF;
    }
}
