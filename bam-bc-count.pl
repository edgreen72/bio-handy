#!/usr/bin/perl

use Bio::DB::Sam;
use Getopt::Std;
use vars qw( $opt_b $opt_q $opt_5 $opt_7 );
use strict;

&init();

my $sam = Bio::DB::Sam->new( -bam   => $opt_b,
			     -expand_flags => 1 );
my ( $feature, $p5_p, $p7_p, $p5, $p7, $comp_bc, $count );
my ( %bc );
my $iterator = $sam->features(-iterator=>1);
#foreach $feature ( $sam->features('match') ) {
while( $feature = $iterator->next_seq ) {
    unless( $feature->unmapped ) {
	if ( $feature->qual < $opt_q ) {
	    next;
	}
	$bc{uc($feature->get_tag_values('BC'))}++;
    }
}

$p5_p = &parse_bc( $opt_5 );
$p7_p = &parse_bc( $opt_7 );

foreach $p5 ( @{ $p5_p } ) {
    foreach $p7 ( @{ $p7_p } ) {
	$comp_bc = "$p5"."$p7";
	if ( $bc{ $comp_bc } ) {
	    $count = $bc{ $comp_bc };
	}
	else {
	    $count = 0;
	}
	printf( "%d %s\n", $count, $comp_bc );
    }
}

sub init {
    getopt( 'b:q:5:7:' );
    unless( -f $opt_b &&
	    -f $opt_5 &&
	    -f $opt_7 ) {
	print( "bam-bc-count.pl -b <bam file> -5 <P5 barcode> -7 <P7 barcodes>\n" );
	print( "                -q <map quality cutoff>\n" );
	print( "Reports a count of the observed number of each barcode combination\n" );
	print( "in the input bam file.\n" );
	exit( 0 );
    }
}

sub parse_bc {
    my $fn = shift;
    my $l;
    my @bc;
    open( BC, $fn ) or die( "$!: $fn\n" );
    while( chomp( $l = <BC> ) ) {
	push( @bc, $l );
    }
    close( BC );
    return \@bc;
}
