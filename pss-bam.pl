#!/usr/bin/perl

use Bio::DB::Sam;
use Getopt::Std;
use vars qw( $opt_f $opt_b $opt_r $opt_l $opt_L $opt_q $opt_F $opt_B $opt_m );
use strict;

my $VERSION = 0.02;

my ( $sam, $feature, $ref, $matches, $query, $i, $rcref, $rcquery, $pair );
my ( @ref, @query, @rref, @rquery );
my (@SUBS, @RSUBS );
my @BASES = ('A', 'C', 'G', 'T' );
my %RC = ('A' => 'T', 'C' => 'G', 'G' => 'C', 'T' => 'A',
	  'a' => 't', 'c' => 'g', 'g' => 'c', 't' => 'a',
	  'N' => 'N', 'n' => 'n', 'X' => 'x', '-' => '-');
	  
&init();

my $sam = Bio::DB::Sam->new( -bam   => $opt_b,
			     -fasta => $opt_f,
			     -expand_flags => 1 );

foreach $feature ( $sam->features('match') ) {
    unless( $feature->unmapped ) {
	($ref, $matches, $query) = $feature->padded_alignment();
	
	### Apply length and map-quality filters
	if ( (length( $query ) < $opt_l) ||
	     (length( $query ) > $opt_L) ||
	     ($feature->qual < $opt_q) ||
	     !&barcode_check($feature->get_tag_values('BC'))  ||
	     (!$opt_m &&  !$feature->get_tag_values('MAP_PAIR')) ) {
	    next;
	}
	
	### It's mapped and passes filters
	### Now, we have to handle paired-end different than merged reads
	if ( $feature->strand == -1 ) {
	    $ref = &revcom( $ref );
	    $query = &revcom( $query );
	}
	@ref   = split( '', $ref );
	@query = split( '', $query );

	### If we have merged reads or if we're looking at the first mate (P5)
	### end
	if ( $opt_m ||
	     $feature->get_tag_values('FIRST_MATE') ) {
	    for( $i = 0; $i < $opt_r; $i++ ) {
		$SUBS[$i]->{$ref[$i]}->{$query[$i]}++;
	    }
	}
	elsif ( $feature->get_tag_values('SECOND_MATE') ) {
	    for( $i = 0; $i < $opt_r; $i++ ) {
		$RSUBS[$i]->{$ref[$i]}->{$query[$i]}++;
	    }
	}   
	
	### Now the other way
#	$rcref   = &revcom( $ref );
#	$rcquery = &revcom( $query );
	####    NOTE: NOT reverse complement. We want to see what's on the top strand
	####          here, too - not the reverse complement of it, which would be
	####          the bottom strand
	@rref = reverse (split( '', $ref ));
	@rquery = reverse (split( '', $query ));

	if ( $opt_m ) {
	    for( $i = 0; $i < $opt_r; $i++ ) {
		$RSUBS[$i]->{$rref[$i]}->{$rquery[$i]}++;
	    }
	}
    }
}



print( "### pss-bam.pl v $VERSION $opt_f\n" );
&output( '### Forward read substitution counts',
	 \@SUBS );

&output( "\n\n### Reverse read substitution counts",
	 \@RSUBS );

0;

# Return true if barcode filter was set and this matches
# Return true if no barcode filter was set
# Return false otherwise
sub barcode_check {
    my $barcode = shift;
    if ( defined( $opt_F ) ) {
	unless( $barcode =~ /^$opt_F/ ) {
	    return 0;
	}
    }
    if ( defined( $opt_B ) ) {
	unless( $barcode =~ /$opt_B$/ ) {
	    return 0;
	}
    }
    return 1;
}

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
    my $l_DEF = 0;
    my $L_DEF = 250000000;
    my $q_DEF = 0;
    getopts( 'f:b:r:l:L:q:F:B:m' );
    unless( -f $opt_f &&
	    -f $opt_b ) {
	print( "pss-bam.pl v $VERSION -f <fasta file> -b <bam file>\n" );
	print( "           -r <region length; default = $r_DEF\n" );
	print( "           -l <min length filter; default = $l_DEF>\n" );
	print( "           -L <max length filter; default = $L_DEF>\n" );
	print( "           -q <map quality filter>\n" );
	print( "           -F <front barcode sequence>\n" );
	print( "           -B <back barcode sequence>\n" );
	print( "           -m <run in merged mode => reads are merged>\n" );
	print( "Uses the Bio::DB::Sam module to make map-damage like data\n" );
	print( "suitable for plotting in gnuplot.\n" );
	exit( 0 );
    }
    unless( defined( $opt_r ) ) {
	$opt_r = $r_DEF;
    }
    unless( defined( $opt_l ) ) {
	$opt_l = $l_DEF;
    }
    unless( defined( $opt_L ) ) {
	$opt_L = $L_DEF;
    }
    unless( defined( $opt_q ) ) {
	$opt_q = $q_DEF;
    }
    unless( defined( $opt_m ) ) {
	$opt_m = 0;
    }
}
