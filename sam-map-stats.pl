#!/usr/local/bin/perl

use strict;
my ( $l,  
     $t, # total lines
     $mapped, # total mapped reads
     $tot_len, # total length of all mapped reads
    );
my ( @sam );

my $id = shift(@ARGV);

unless ( defined( $id ) ) {
    $id = "NULL";
}

while( chomp( $l = <STDIN> ) ) {
    $t++;
    @sam = split( "\t", $l );
    if ( $sam[2] ne "*" ) { # It mapped
	$mapped++;
	$tot_len += length($sam[9]);
    }
}

printf( "%s %d %.3f %.1f\n", $id, $t, ($mapped/$t), ($tot_len/$mapped));

