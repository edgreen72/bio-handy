#!/usr/local/bin/perl

use Getopt::Std;
use vars qw( $opt_d $opt_s );
use strict;

&init();
my ( $inx2id_p, $inx2count_p, $inx, $id );
$inx2id_p    = &parse_SampleSheet( $opt_s );
$inx2count_p = &parse_demultSummary( $opt_d );

foreach $inx ( sort { $inx2count_p->{$a} <=> $inx2count_p->{$b} } 
	       keys %{ $inx2count_p } ) {
    if ( $inx2id_p->{$inx} ) {
	$id = $inx2id_p->{$inx};
    }
    else {
	$id = "NA";
    }
    printf( "%s %s %d\n", $id, $inx, $inx2count_p->{$inx} );
}

sub init {
    getopts( 'd:s:' );
    unless( -f $opt_d &&
	    -f $opt_s ) {
	print( "MiSeqRunCheck.pl -s <SampleSheet.csv> -d <DemultiplexSummaryF1L1.txt\n" );
	print( "Reports the IDs and counts of the most common indeces seen in a run\n" );
	print( "and writes NA if the index was not listed in the SampleSheet.csv\n" );
	exit( 0 );
    }
}

sub parse_demultSummary {
    my $fn = shift;
    my ( $l, $inx, $rc_inx, $count );
    my ( %inx2count );

    open( DS, $fn ) or die( "$!: $fn\n" );
    chomp( $l = <DS> );
    until ( $l =~ /^Index/ ) {
	chomp( $l = <DS> );
    }
    while( chomp( $l = <DS> ) ) {
	if ( $l =~ /^Index2/ ) {
	    close( DS );
	    return \%inx2count;
	}
	($inx, $rc_inx, $count) = split( ' ', $l );
	$inx2count{ $inx } = $count;
    }
    close( DS );
    return \%inx2count;
}
    
    
sub parse_SampleSheet {
    my $fn = shift;
    my ( $l, $i, $id_inx, $inx_inx );
    my ( @cats, @sample );
    my ( %inx2id );

    open( SS, $fn ) or die( "$!: $fn\n" );
    chomp( $l = <SS> );
    until ( $l =~ /^\[Data\]/ ) {
	chomp( $l = <SS> );
    }
    chomp( $l = <SS> );
    @cats = split( ',', $l );
    for( $i = 0; $i <= $#cats; $i++ ) {
	if ( $cats[$i] eq 'Sample_ID' ) {
	    $id_inx = $i;
	}
	elsif ( $cats[$i] eq 'index' ) {
	    $inx_inx = $i;
	}
    }

    while( chomp( $l = <SS> ) ) {
	@sample = split( ',', $l );
	$inx2id{ $sample[ $inx_inx ] } = $sample[ $id_inx ];
    }
    close( SS );
    return \%inx2id;
}
	
