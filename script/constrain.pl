#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use Bio::Phylo::IO 'parse_tree';

# process command line arguments
my ( $fasta, $table );
GetOptions(
	'fasta=s'  => \$fasta,
	'table=s'  => \$table,
);

my %exemplars;
{
	open my $fh, '<', $table or die $!;
	while(<$fh>) {
		chomp;
		my ( $tree, $species ) = split /\t/, $_;
		$species =~ s/_/ /g;
		$exemplars{$species} = $tree;	
	}
}

my %pair;
{
	open my $fh, '<', $fasta or die $!;
	while(<$fh>) {
		chomp;
		if ( /^>(.+)/ ) {
			my $defline = $1;
			my @parts   = split /\|/, $defline;
			my $id      = $parts[0];
			my $species = $parts[1];
			my $tree    = $exemplars{$species};
			$id = substr $id, 0, 10;			
			$pair{$tree} = [] if not $pair{$tree};
			push @{ $pair{$tree} }, $id;
		}
	}
}

print '(';
print join ',', map { '(' . join(',',@$_) . ')' } values %pair;
print ');';