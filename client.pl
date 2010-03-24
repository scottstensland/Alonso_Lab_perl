
use Alonso_Lab::Motif_Discovery;
use Data::Dumper;
use strict;
use warnings;

package main;  #  just use default namespace here in the client


my $motif_discovery = new Motif_Discovery();

#print Dumper ( $item );


my %hash_stuff;

$hash_stuff{'apple'} = 5;
$hash_stuff{'pear'} = 6;

$motif_discovery->firstName('Scott', \%hash_stuff );


print Dumper( \%hash_stuff);


$motif_discovery->lastName('Stensland');

$motif_discovery->print();