#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use Bio::Phylo::IO 'parse_tree';
use Bio::Phylo::Util::Logger ':levels';

# process command line arguments
my $verbosity = WARN;
my $intree;
GetOptions(
	'intree=s' => \$intree,
	'verbose+' => \$verbosity,
);

# instantiate logger
my $log = Bio::Phylo::Util::Logger->new(
	'-level' => $verbosity,
	'-class' => 'main',
);

# read the tree
$log->info("going to read tree from $intree");
my $tree = parse_tree(
	'-format' => 'newick',
	'-file'   => $intree,
);
$log->info("instantiated tree $tree");

# every tree has a single outgroup (I think!), so
# we first want to traverse to the ingroup, and from
# there get the leftmost and rightmost terminal
my @children = @{ $tree->get_root->get_children };
if ( @children == 2 ) {
	my ( $outgroup, $mrca ) = sort { scalar(@{$a->get_terminals}) <=> scalar(@{$b->get_terminals}) } @children;
	print $intree, "\t", $mrca->get_leftmost_terminal->get_name, "\n";
	print $intree, "\t", $mrca->get_rightmost_terminal->get_name, "\n";
}
else {
	$log->error("unexpected number of root children (".scalar(@children).") in $intree");
	die;
}