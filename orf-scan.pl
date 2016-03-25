#!/usr/bin/perl

### (c) 2016 - UC Regents
### Author: Ed Green
###         Dept. of Biomolecular Engineering
###         UC Santa Cruz

use Bio::SeqIO;
use Getopt::Std;
use vars qw( $opt_f $opt_s $opt_t $opt_m );                                                                                   
my $VERSION = 1.0;

### Process command-line arguments, set up Bio::SeqIO object
### start and stop codons
my ($bio_io, $starts_p, $terms_p) = &init();                                                                                  

my ( $seq_o, $next_start, $next_term );                                                                                       
my ( @starts, @terms );                                                                                                       
my ( %orfs );                                                                                                                 


### Go through each input sequence in the input fasta database
while( $seq_o = $bio_io->next_seq() ) {
    
    ### For each reading frame 0, 1, & 2
    for($i = 0; $i<= 2; $i++ ) {                                                                                              

	### Find 0-index positions of all start and stop codons
	### in this sequence, in this frame, and return sorted list
        @starts = &find_codons( $seq_o, $i, $starts_p );                                                                      
        @terms  = &find_codons( $seq_o, $i, $terms_p );                                                                       

	### Initialize next_start and next_term
        $next_start = shift( @starts );
        $next_term  = -1;

	### Find the next termination codon, add it to %orfs and then
	### skip to the next start codon past that termination codon
	### continue until there are no more start codons
        while( ($#starts >= 0) &&
               ($#terms  >= 0) ) {
            while ( ($#terms >= 0) &&                                                                     
                    ($next_term <= $next_start) ) {
                $next_term = shift(@terms);
            }
            if ( $next_start < $next_term ) {
                if ( ($next_term - $next_start) >= $opt_m ) {
                    push( @{$orfs{'Forward'}->{$i}}, [$next_start, $next_term] );
                }
                while( ($#starts >= 0) &&
                       ($next_start <= $next_term) ) {
                    $next_start = shift( @starts );
                }
            }
        }
    }

    ### same as above, but with reverse complement of sequence
    $rev_seq_o = $seq_o->revcom();
    $length = $rev_seq_o->length();
    for($i = 0; $i<= 2; $i++ ) {
        @starts = &find_codons( $rev_seq_o, $i, $starts_p );
        @terms  = &find_codons( $rev_seq_o, $i, $terms_p );

        $next_start = shift( @starts );
        $next_term  = -1;
        while( ($#starts >= 0) &&
               ($#terms  >= 0) ) {
            while ( ($#terms >= 0) &&
                    ($next_term <= $next_start) ) {
                $next_term = shift(@terms);
            }
            if ( $next_start < $next_term ) {
                if ( ($next_term - $next_start) >= $opt_m ) {
                    push( @{ $orfs{'Reverse'}->{$i}}, 
                          [$next_start, $next_term] );
                }
                while( ($#starts >= 0) &&
                       ($next_start <= $next_term) ) {
                    $next_start = shift( @starts );
                }
            }

        }
    }

    &make_output( \%orfs, $seq_o );
    %orfs = ();
}

sub init {
    my $bio_io;
    my @starts;
    my @terms;
    my $s_DEF = 'ATG';
    my $t_DEF = 'TAG:TAA:TGA';
    my $m_DEF = 27;
    getopts( 'f:s:t:m:' );
    unless( -f $opt_f ) {
        print( "orf-scan.pl Version $VERSION -f <fasta genome> -s <STARTs> -t <TERMs> -m <MIN>\n" );
        print( "Scans the fasta genome sequence for the input start and termination\n" );
        print( "codons and outputs all observed ORFs greater than MIN length.\n" );
        print( "STARTs and TERMs should be a colon delimited list. For example:\n" );
        print( "-s $s_DEF -t $t_DEF\n" );
        print( "Writes an output file with the folowing columns\n" );
        print( "1. ID\n" );
        print( "2. Strand (Forward or Reverse)\n" );
        print( "3. Frame (0, 1, or 2)\n" );
        print( "4. Position of first base of start codon\n" );
        print( "5. Position of first base of stop codon\n" );
        print( "-s DEFAULT = $s_DEF\n" );
        print( "-t DEFAULT = $t_DEF\n" );
        print( "-m DEFAULT = $m_DEF\n" );
        exit( 0 );
    }
    unless( $opt_s ) {
        $opt_s = $s_DEF;
    }
    unless( $opt_t ) {
        $opt_t = $t_DEF;
    }
    unless( $opt_m ) {
        $opt_m = $m_DEF;
    }
    my $bio_io = Bio::SeqIO->new( '-file' => $opt_f, '-format' => 'fasta' );
    @starts = split( ':', $opt_s );
    @terms  = split( ':', $opt_t );

    return ( $bio_io, \@starts, \@terms );
}
    
### make_output
### This is called on each sequence
sub make_output {
    my $orfs_p = shift;
    my $seq_o  = shift;
    my $rc_seq_o = $seq_o->revcom();
    my $id = $seq_o->display_id();
    my $length = $seq_o->length();
    my @in_orf = 0 x $length;
    my ( $strand, $frame, $orf_p, $translation, $i, $total );
    foreach $strand ( keys %{ $orfs_p } ) {
        foreach $frame ( keys %{ $orfs_p->{$strand}} ) {
            foreach $orf_p (@{ $orfs_p->{$strand}->{$frame} }) {
                if ( $strand eq 'Forward' ) {
                    $translation = $seq_o->trunc($orf_p->[0]+1, $orf_p->[1]+3)->
                        translate()->seq;
                    printf( "%s\t%s\t%d\t%d\t%d\t%s\n",
                            $id, $strand, $frame, $orf_p->[0]+1, $orf_p->[1]+3, $translation );
                    for ( $i = $orf_p->[0]; $i <= $orf_p->[1]; $i++ ) {
                        $in_orf[$i] = 1;
                    }
                }
                else {
                    $translation = $rc_seq_o->trunc($orf_p->[0]+1, $orf_p->[1]+3)->
                        translate()->seq;
                    printf( "%s\t%s\t%d\t%d\t%d\t%s\n",
                            $id, $strand, $frame, 
                            $length-$orf_p->[0], $length-$orf_p->[1]+2, $translation );
                    for ( $i = $orf_p->[0]; $i <= $orf_p->[1]; $i++ ) {
                        $in_orf[$i] = 1;
                    }
                }
            }
        }
    }
    for( $i = 0; $i < $length; $i++ ) {
        $total += $in_orf[$i];
    }
    printf( "# %d of %d in ORFs in %s\n", $total, $length, $id );
}

sub find_codons {
    my $seq_o = shift;
    my $i     = shift;
    my $codons_p = shift;
    my $seq_str = $seq_o->seq;
    my ( $codon, $pos );
    my ( @all_pos, @pos );
    foreach $codon ( @{ $codons_p } ) {
        $pos = index( $seq_str, $codon, $i );
        while( $pos > 0 ) {
            push( @all_pos, $pos );
            $pos = index( $seq_str, $codon, $pos+1 );
        }
    }
    foreach $pos ( @all_pos ) {
        if ( $pos%3 == $i ) {
            push( @pos, $pos );
        }
    }
    return sort {$a <=> $b} @pos;
}
 
