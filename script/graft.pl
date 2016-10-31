#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use Bio::Phylo::IO 'parse_tree';
use Bio::Phylo::Util::Logger ':levels';

# process command line arguments
my $verbosity = WARN;
my ( $backbone, $dir, $table );
GetOptions(
	'verbose+'   => \$verbosity,
	'backbone=s' => \$backbone,
	'dir=s'      => \$dir,
	'table=s'    => \$table,
);

# instantiate helper objects
my $log = Bio::Phylo::Util::Logger->new(
	'-level' => $verbosity,
	'-class' => 'main',
);

$log->info("going to read backbone $backbone");
my $bb = parse_tree(
	'-format' => 'newick',
	'-file'   => $backbone,
);

my %clade;
{
	$log->info("going to read mapping $table");
	open my $fh, '<', $table or die $!;
	while(<$fh>) {
		chomp;
		my ( $clade, $species ) = split /\t/, $_;
		$clade{$clade} = [] if not $clade{$clade};
		push @{ $clade{$clade} }, $species;
	}
}

# start grafting clades
for my $clade ( sort { $a cmp $b } keys %clade ) {
	
	# read the clade tree, fetch the ingroup
	$log->info("going to graft clade $clade");
	my $tree = parse_tree(
		'-format' => 'newick',
		'-file'   => "${dir}/${clade}",
	);
	my ( $ingroup ) = sort { scalar(@{$b->get_terminals}) <=> scalar(@{$a->get_terminals}) } @{ $tree->get_root->get_children };

	# fetch the equivalent exemplar clade
	my @tips = map { $bb->get_by_name($_) } @{ $clade{$clade} };
	$log->info("found ".scalar(@tips)." exemplar tips");
	my $mrca = $bb->get_mrca(\@tips);
	
	# do the graft
	$ingroup->set_branch_length( $mrca->get_branch_length );
	$ingroup->set_parent( $mrca->get_parent );
	$ingroup->visit_depth_first( '-pre' => sub { $bb->insert(shift) } );
	$bb->prune_tips( $mrca->get_terminals );
}

print $bb->to_newick;