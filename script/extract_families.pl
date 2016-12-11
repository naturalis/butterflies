#!/usr/bin/perl
use strict;
use Getopt::Long;
use Bio::Phylo::IO 'parse_tree';

# process command line arguments
my $dir;
my $prefix    = 'RAxML_bipartitions.';
my $extension = 'nwk';
GetOptions(
	'dir=s'       => \$dir,
	'prefix=s'    => \$prefix,
	'extension=s' => \$extension,
);

# start reading directories
opendir my $dh, $dir or die $!;
while( my $entry = readdir $dh ) {
	if ( $entry =~ /$prefix(.+)\.$extension/ ) {
		my $family = $1;
		$family =~ s/_rooted$//g;
		if ( $family =~ /[A-Z][a-z]+/ ) {

			warn "reading $dir/$entry for clade $family";
			my $ft = parse_tree(
				'-format' => 'newick',
				'-file'   => "$dir/$entry",
			);
			my ($ingroup) = sort {
				scalar @{ $b->get_terminals }
				<=>
				scalar @{ $a->get_terminals }
			} @{ $ft->get_root->get_children };
			for my $tip ( @{ $ingroup->get_terminals } ) {
				my $name = $tip->get_name;
				print $name, "\t", $family, "\n";
			}
		}
	}
}