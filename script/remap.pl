#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use Bio::Phylo::IO 'parse_tree';

# process command line arguments
my ( $fasta, $intree );
GetOptions(
	'fasta=s'  => \$fasta,
	'intree=s' => \$intree,
);

my $tree = parse_tree(
	'-format' => 'newick',
	'-file'   => $intree,
);

my %taxon;
open my $fh, '<', $fasta or die $!;
while(<$fh>) {
	chomp;
	if ( /^>(.+)/ ) {
		my $defline = $1;
		my @parts = split /\|/, $defline;
		my $id = $parts[0];
		my $taxon = $parts[1];
		$id = substr $id, 0, 10;
		$taxon{$id} = $taxon;
	}
}

for my $tip ( @{ $tree->get_terminals } ) {
	my $id = $tip->get_name;
	if ( $taxon{$id} ) {
		$tip->set_name( $taxon{$id} );
	}
	else {
		warn $id;
	}
}

print $tree->to_newick;