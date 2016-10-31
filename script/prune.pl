#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use Bio::Phylo::Factory;
use Bio::Phylo::IO qw'parse_tree unparse';
use Bio::Phylo::Util::Logger ':levels';

# process command line arguments
my $verbosity = WARN;
my ( $spreadsheet, $intree, $missing );
GetOptions(
	'spreadsheet=s' => \$spreadsheet,
	'intree=s'      => \$intree,
	'verbose+'      => \$verbosity,
	'missing'       => \$missing,
);

# instantiate logger
my $log = Bio::Phylo::Util::Logger->new(
	'-level' => $verbosity,
	'-class' => 'main',
);

# read tree
$log->info("going to read tree $intree");
my $tree = parse_tree(
	'-format' => 'newick',
	'-file'   => $intree,
);

# collect tips to keep
my @keep;
my @missing;
$log->info("going to keep taxa from $spreadsheet");
open my $fh, '<', $spreadsheet or die $!;
while(<$fh>) {
	chomp;
	my @line = split /\t/, $_;
	if ( $line[1] ne 'Species' ) {
		my $name = $line[1];
		$name =~ s/ /_/g;
		my $tip = $tree->get_by_name($name);
		if ( $tip ) {
			push @keep, $tip;
		}
		else {
			$log->warn("$name is not in tree");
			push @missing, $name;	
		}
	}
}

$log->info("found ".scalar(@keep)." tips to keep");
$log->info(scalar(@missing)." species were absent from tree");

# write out the missing 
if ( $missing ) {
	print join "\n", @missing;
}
else {
	$tree->keep_tips(\@keep);
	
	# convert to nexus
	my $fac = Bio::Phylo::Factory->new;
	my $proj = $fac->create_project;
	my $forest = $fac->create_forest;
	$forest->insert($tree);
	my $taxa = $forest->make_taxa;
	$proj->insert($taxa);
	$proj->insert($forest);
	print $proj->to_nexus;
}