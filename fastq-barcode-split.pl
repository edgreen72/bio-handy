#!/usr/local/bin/perl

use Getopt::Std;
use vars qw( $opt_f $opt_r $opt_l $opt_o );
use strict;

my $TOTAL = 0;
my $READ_PAIRS;
my ( $for_fq, $rev_fq, $for_bc, $rev_bc );
my ( %BC_counts );
my %RC = ( 'A' => 'T', 
	   'C' => 'G',
	   'G' => 'C',
	   'T' => 'A',
	   'N' => 'N',
	   'X' => 'X',
	   'a' => 't',
	   'c' => 'g',
	   'g' => 'c',
	   't' => 'a',
	   'n' => 'n',
	   'x' => 'x' );
my %BC2FH; # barcode 2 filehandle hash

&init();

if ( $opt_f =~ /\.gz$/ ) {
    open( FOR, "zcat $opt_f|" ) or die( "$!: $opt_f\n" );
}
else {
    open( FOR, $opt_f ) or die( "$!: $opt_f\n" );
}

### Did user specify reverse reads file?
if ( -f $opt_r ) {

    $READ_PAIRS = 1;
    if ( $opt_r =~ /\.gz$/ ) {
	open( REV, "zcat $opt_r|" ) or die( "$!: $opt_r\n" );
    }
    else {
	open( REV, $opt_r ) or die( "$!: $opt_r\n" );
    }
}
else {
    $READ_PAIRS = 0;
}

### Go through the input, sorting by barcodes
if ( $READ_PAIRS ) {
    while( ($for_fq = &read_next_fq( \*FOR )) &&
	   ($rev_fq = &read_next_fq( \*REV )) ) {
	$for_bc = &read_barcode( $for_fq, $opt_l );
	$rev_bc = &read_barcode( $rev_fq, $opt_l );
	
	$BC_counts{ $for_bc }{$rev_bc }++;
	$TOTAL++;
    }

    close( FOR );
    close( REV );
}

else {
    ### No pairs, so this must be a single "merged" file.
    ### Get the barcodes from the front and the revcom of the end
    while( $for_fq = &read_next_fq( \*FOR )) {
	($for_bc, $rev_bc) = &read_barcodes( $for_fq, $opt_l );
	$BC_counts{ $for_bc }{ $rev_bc }++;
	$TOTAL++;
	if ( $opt_o ) {
	    &write_bc_trunc_fq( $opt_l, $for_fq, $for_bc, $rev_bc );
	}
    }
    close( FOR );
}

foreach $for_bc ( keys %BC_counts ) {
    foreach $rev_bc ( keys %{ $BC_counts{$for_bc} } ) {
	printf( "%s %s %d %0.4f\n", 
		$for_bc, 
		$rev_bc, 
		$BC_counts{$for_bc}{$rev_bc},
		$BC_counts{$for_bc}{$rev_bc}/$TOTAL,
	    );
    }
}

sub read_barcode {
    my $fq_p        = shift;
    my $barcode_len = shift;
    my $barcode;
    $barcode = substr( $fq_p->{'seq'}, 0, $barcode_len );
    return $barcode;
}

sub read_barcodes {
    my $fq_p        = shift;
    my $barcode_len = shift;
    my $f_barcode;
    my $r_barcode;
    $fq_p->{'rc_seq'} = &revcom( $fq_p->{'seq'} );
    $f_barcode = substr( $fq_p->{'seq'}, 0, $barcode_len );
    $r_barcode = substr( $fq_p->{'rc_seq'}, 0, $barcode_len );
    return ($f_barcode, $r_barcode);
}

sub revcom {
    my $seq = shift;
    my @seq = split('', $seq);
    my $rcseq;
    
    $rcseq = join( '', (map { &rc_base($_) } reverse( @seq )) );
    return $rcseq;
}

sub rc_base {
    my $b = shift;
    if ( defined( $RC{$b} ) ) {
	return $RC{$b};
    }
    return $b;
}

sub read_next_fq {
    my $fh = shift;
    my %FQ;
    my $l;
    if ( chomp( $l = <$fh> ) ) { # got header
	if ( $l =~ /^\@/ ) { # looks like header
	    $FQ{ 'header' } = $l;
	    if ( chomp( $l = <$fh> ) ) { # got sequence
		$FQ{'seq'} = $l;
		if ( chomp( $l = <$fh> ) ) { # got qual header 
		    $FQ{'qheader'} = $l;
		    if ( chomp( $l = <$fh> ) ) { # got qual
			$FQ{'qual'} = $l;
			return \%FQ;
		    }
		}
	    }
	}
    }
    return 0; # something didn't work
}

sub write_bc_trunc_fq {
    my $trunc_len = shift;
    my $fq_p      = shift;
    my $for_bc    = shift;
    my $rev_bc    = shift;

    my $comp_bc = $for_bc.$rev_bc;
    my $fh;
    my ( $trunc_seq, $trunc_qual );

    ### Get filhandle, make if necessary
    unless ( defined( $BC2FH{ $comp_bc } ) ) {
	open( $BC2FH{ $comp_bc }, ">", $opt_o.'.'.$comp_bc.'.fq' );
    }
    $fh = $BC2FH{$comp_bc};

    ### Truncate the sequence by the barcode length at the beginning
    $trunc_seq  = substr( $fq_p->{'seq'}, $trunc_len );
    $trunc_qual = substr( $fq_p->{'qual'}, $trunc_len );
    $fq_p->{'seq'}  = $trunc_seq;
    $fq_p->{'qual'} = $trunc_qual;

    ### If we have only merged (not pairs), then truncate the end, too
    unless( $READ_PAIRS ) {
	$trunc_seq  = substr( $fq_p->{'seq'},  0, length($fq_p->{'seq'})  - $trunc_len );
	$trunc_qual = substr( $fq_p->{'qual'}, 0, length($fq_p->{'qual'}) - $trunc_len );
	$fq_p->{'seq'}  = $trunc_seq;
	$fq_p->{'qual'} = $trunc_qual
    }

    print $fh ($fq_p->{'header'}, "\n",
	       $fq_p->{'seq'}, "\n",
	       $fq_p->{'qheader'}, "\n",
	       $fq_p->{'qual'}, "\n" );
}
    
    

sub init {
    my $l_DEF = 5;
    getopts( 'f:r:l:o:' );
    unless( -f $opt_f ) {
	print( "fastq-barcode-split.pl -f <forward fastq file> -r <reverse fastq file> -l [BARCODE length]\n" );
	print( "                       -o <root name for output file(s) to make>\n" );
	print( "The program reads one or a pair of fastq files specified by the -f and -r options.\n" );
	print( "If no -r file is given, it's ass-u-med that the single input file has reads with\n" );
	print( "no adapters and barcodes at the very beginning and end of the reads.\n" );
	print( "If both -f and -r files are given, they are ass-u-med to be the forward and reverse\n" );
	print( "reads from a paired-end run. Reads must be in the same order.\n" );
	print( "It takes the first -l bases of the forward and reverse reads and treats them as a\n" );
	print( "composite barcode. It reports the number of read pairs with each observed barcode.\n" );
	print( "For example, a read pair like this, using -l 5\n" );
	print( '@M00160:20:000000000-ANYJ6:1:1101:10222:1120 1:N:0:8\n' );
	print( "GATGAGGCTTAGAAGCAAGATCGGAAGAGCACACGTCTGAACTCCAGTCACTTGGTCAACATCGTTTGCCGTCTTC\n" );
	print( "+\n" );
	print( 'CC@CCGGGEFFACFFFGC6,,CCCCC+6CE@<FFDCFFG9C@F9<@CFDAFCG<,;CAF,,;@,C,C,EF+CEF7C\n' );
	print( "[other read]\n" );
	print( '@M00160:20:000000000-ANYJ6:1:1101:10222:1120 2:N:0:8\n' );
	print( "TGCTTCTAAGCCTAATCAGATCGGAAGAGCGTCGTGTAGNGAAAGANTNNAGNTNTCGNTNGGNGCCGTATCATTA\n" );
	print( "+\n" );
	print( 'CCCCCGGGGG?8F,;CFFEDFF,@FGD8@D@F,66BCEC#6,CC@<#:##,,#6#::C#,#:,#:C,@C@,,<,,<', "\n" );
	print( "Would be classified as GATGA.TGCTT\n" );
	exit( 0 );
    }
    unless( defined( $opt_l ) ) {
	$opt_l = $l_DEF;
    }
}
