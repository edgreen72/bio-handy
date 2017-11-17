#!/usr/local/bin/perl

use Bio::AlignIO;
use Bio::SeqIO;
use Getopt::Std;
use vars qw( $opt_a $opt_c $opt_I );
use strict;

my $aln_o = &init();

my $consensus = $aln_o->consensus_string( $opt_c );

$consensus =~ s/\?/N/g;

my $seq_io = Bio::SeqIO->new( '-format' => 'fasta' );

my $seq_o = Bio::Seq->new( '-seq' => $consensus,
			   '-display_id' => $opt_I );

$seq_io->write_seq( $seq_o );

sub init {
    my $c_DEF = 50;
    my $I_DEF = 'Consensus';
    getopts( 'a:c:I:' );
    unless( -f $opt_a ) {
	print( "align2consensus.pl\n" );
	print( "  -a <alignment file; clustalw format>\n" );
	print( "  -c <percent ID for consensus call; default $c_DEF>\n ");
	print( "  -I <Identifier to use; default $I_DEF>\n" );
	exit( 0 );
    }

    return Bio::AlignIO->new( '-file' => $opt_a,
			      '-format' => 'clustalw')->next_aln;
}

    
