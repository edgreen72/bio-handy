#!/usr/bin/perl

use Bio::DB::Sam;
use Getopt::Std;
use vars qw( $opt_f $opt_b $opt_r $opt_l $opt_L $opt_q $opt_F $opt_B $opt_m );
use strict;

my $VERSION = 0.05;
my $DEBUG = 1;
my $CONTEXT = 2;
my ( $sam, $feature, $ref, $matches, $query, $i, $rcref, $rcquery, $pair, $length );
my ( $ref_id, $ref_start, $ref_end, $up_dna, $down_dna, $tmp_dna );
my ( @ref, @query, @rref, @rquery );
my (@SUBS, @RSUBS );
my @BASES = ('A', 'C', 'G', 'T' );
my $up_context_p;  # ->[-CONTEXT] -> {'A' => counts, 'C' => counts, ...
                   #   [-CONTEXT + 1] -> {'A' => counts, 'C' => counts, ...
my $down_context_p;# ->[1] -> {'A' => counts, 'C' => counts, ...
                   #   [2] -> {'A' => counts, 'C' => counts, ...
                   #...[CONTEXT] -> {'A' => counts, 'C' => counts, ...
my %RC = ('A' => 'T', 'C' => 'G', 'G' => 'C', 'T' => 'A',
	  'a' => 't', 'c' => 'g', 'g' => 'c', 't' => 'a',
	  'N' => 'N', 'n' => 'n', 'X' => 'x', '-' => '-');
	  
&init();

my $sam = Bio::DB::Sam->new( -bam   => $opt_b,
			     -fasta => $opt_f,
			     -expand_flags => 1 );

#foreach $feature ( $sam->features('match') ) {
my $iterator = $sam->features(-iterator=>1);
while ( $feature = $iterator->next_seq ) {
    unless( $feature->unmapped ) {
	($ref, $matches, $query) = $feature->padded_alignment();
	if ( $opt_m ) {
	    $length = length($query);
	}
	else {
	    $length = abs($feature->isize); # size is reported as negative sometimes
	}
	### Apply filters
	if ( ($length < $opt_l) || # too small
	     ($length > $opt_L) || # too big
	     ($feature->qual < $opt_q) || # too low map quality
	     !&barcode_check($feature->get_tag_values('BC'))  || # wrong barcode
	     (!$opt_m &&  !$feature->get_tag_values('PAIRED')) ) {
	    next;
	}
	
	### It's mapped and passes filters
	$ref_id    = $feature->seq_id;
	$ref_start = $feature->start;
	$ref_end   = $feature->end;
	$up_dna   = uc($sam->seq($ref_id, ($ref_start - $CONTEXT), ($ref_start - 1)));
	$down_dna = uc($sam->seq($ref_id, ($ref_end + 1), ($ref_end + $CONTEXT)));
	unless( length($up_dna) == $CONTEXT ) {
	    $up_dna = 'N' x $CONTEXT;
	}
	unless( length($down_dna) == $CONTEXT ) {
	    $down_dna = 'N' x $CONTEXT;
	}

	### If this is a single, merged sequence or if it's the forward (aka P5
	### aka FIRST_MATE), then revcom if it's aligned to the minus strand
	### so that we get the actual sequence as it was read against the
	### minus strand genome sequence.
	if ( $opt_m ||
	     ($feature->get_tag_values('FIRST_MATE')) ) {
	    if ( $feature->strand == -1 ) {
		$ref   = &revcom( $ref );
		$query = &revcom( $query );
		$tmp_dna  = $up_dna;
		$up_dna   = &revcom( $down_dna );
		$down_dna = &revcom( $up_dna );
	    }
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
	    &add_up_context( $up_dna );
	}

	elsif ( $feature->get_tag_values('SECOND_MATE') ) {

	    if ( $feature->get_tag_values('REVERSED') !=
		 $feature->get_tag_values('M_REVERSED') ) { # Innie's only

		if ( $feature->get_tag_values('REVERSED') ) {
		    # P7 read is revcom in its alignment to the reference.
		    # Keep it this way, since that is the strand of the P5
		    # alignment. But, reverse it so we have the sequence
		    # from the end of the read toward the interior.
		    @rref   = reverse (split( '', $ref ));
		    @rquery = reverse (split( '', $query ));
		    for( $i = 0; $i < $opt_r; $i++ ) {
			$RSUBS[$i]->{$rref[$i]}->{$rquery[$i]}++;
		    }
		    &add_down_context( $down_dna );
		}
		else { # P7 read is NOT revcom in its alignment.
		       # We want to be on the same strand as the P5
		       # read, which must have been revcom. So, 
		       # revcom, then reverse to get the 3' to 5'
		       # (from end into interior) of the sequence
		       # on the same strand as the P5 read
		    $ref   = &revcom( $ref );
		    $query = &revcom( $query );
		    @rref   = reverse split( '', $ref );
		    @rquery = reverse split( '', $query );
		    for( $i = 0; $i < $opt_r; $i++ ) {
			$RSUBS[$i]->{$rref[$i]}->{$rquery[$i]}++;
		    }
		    &add_down_context( &revcom($up_dna) );
		}
	    }
	}
	elsif ( $opt_m ) {
	    @rref = reverse (split( '', $ref ));
	    @rquery = reverse (split( '', $query ));
	    for( $i = 0; $i < $opt_r; $i++ ) {
		$RSUBS[$i]->{$rref[$i]}->{$rquery[$i]}++;
	    }
	    if ( $feature->strand == -1 ) {
		&add_down_context( &revcom($up_dna) );
	    }
	    else {
		&add_down_context( $down_dna );
	    }
	}
    }
}



print( "### pss-bam.pl v $VERSION\n" );
print( "### $opt_f\n" );
print( "### $opt_b\n" );
&output_beginning( '### Forward read substitution counts and base context',
		   \@SUBS, $up_context_p );

&output_end(       "\n\n### Reverse read substitution counts and base context",
		   \@RSUBS, $down_context_p );
		 

0;

# Takes a string of DNA sequence, $CONTEXT nt long
# Increments bases in $up_context_p
# keys of $up_context_p are negative numbers (upstream)
sub add_up_context {
    my $up_dna       = shift;
    my $i;
    my @up_dna = split( '', $up_dna );
    for( $i = -$CONTEXT; $i < 0; $i++ ) {
	$up_context_p->{ $i }->{ $up_dna[ $i + $CONTEXT ] }++;
    }
}

# Takes a string of DNA sequence, $CONTEXT nt long
# Increments bases in $down_context_p
# keys of $down_context_p are positive numbers (downstream)
sub add_down_context {
    my $down_dna = shift;
    my $i;
    my @down_dna = split( '', $down_dna );
    for( $i = 0; $i < $CONTEXT; $i++ ) {
	$down_context_p->{$i+1}->{ $down_dna[ $i ] }++;
    }
}

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

sub output_beginning {
    my $header_string = shift;
    my $s_p = shift; # substitions
    my $c_p = shift; # base context
    my ( $i, $ref_base, $query_base, $count );
    print( "$header_string\n" );

    for( $i = -$CONTEXT; $i < 0; $i++ ) {
	printf( "%d %d 0 0 0 0 %d 0 0 0 0 %d 0 0 0 0 %d\n", 
		$i, 
		$c_p->{$i}->{'A'},
		$c_p->{$i}->{'C'},
		$c_p->{$i}->{'G'},
		$c_p->{$i}->{'T'} );
    }

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
		printf(" %d", $count);
	    }
	}
	print( "\n" );
    }
}

sub output_end {
    my $header_string = shift;
    my $s_p = shift; # substitions
    my $c_p = shift; # base context
    my ( $i, $ref_base, $query_base, $count );
    print( "$header_string\n" );

#    for( $i = 0; $i < $opt_r; $i++ ) {
    for( $i = $opt_r - 1; $i >= 0; $i-- ) {
	print( "$i" );
	foreach $ref_base (@BASES) {
	    foreach $query_base (@BASES) {
		if ( $s_p->[$i]->{$ref_base}->{$query_base} ) {
		    $count = $s_p->[$i]->{$ref_base}->{$query_base};
		}
		else {
		    $count = 0;
		}
		printf(" %d", $count);
	    }
	}
	print( "\n" );
    }
    for( $i = 1; $i <= $CONTEXT; $i++ ) {
	printf( "%d %d 0 0 0 0 %d 0 0 0 0 %d 0 0 0 0 %d\n", 
		$i, 
		$c_p->{$i}->{'A'},
		$c_p->{$i}->{'C'},
		$c_p->{$i}->{'G'},
		$c_p->{$i}->{'T'} );
    }


}

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
