#!/usr/bin/perl

#initialize
$l = <STDIN>;
@el = split( "\t", $l );
$longest_orf_len = $el[4] - $el[3];
$longest_orf     = $l;
$last_id         = $el[0];

while( chomp( $l = <STDIN> ) ) {
    @el = split( "\t", $l );
    if ( $el[0] eq $last_id ) {
	$len = $el[4] - $el[3];
	if ( $len > $longest_orf_len ) {
	    $longest_orf_len = $len;
	    $longest_orf      = $l;
	}
    }
    else {
	print( "$longest_orf\n" );
	$longest_orf_len          = $el[4] - $el[3];
	$longest_orf              = $l;
	$last_id                  = $el[0];
    }
}

print ("$longest_orf\n");

