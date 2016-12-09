#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use Bio::Phylo::Factory;
use Bio::Phylo::IO 'parse_tree';
use Bio::Phylo::Util::Logger ':levels';

# The rationale is as follows: 
# - I have downloaded a newick file from TimeTree that gives the estimates of node
#   ages for the families in the Lepidoptera. Some of these are primary, others 
#   secondary calibration points. In any case, for every tip (which is a family), its
#   branch length is the STEM age of that family.
# - For every family, I will check if we have that family also as a clade tree from 
#   Mutanen et al. If we do, I will get the leftmost and rightmost tips from the ingroup
#   of the clade tree, find the equivalent ones in the grafted trees, get their MRCA and
#   then apply the age constraint TO THE PARENT NODE OF THE MRCA (it being a stem age
#   estimate, and the MRCA is the crown node).


# process command line arguments
my $verbosity = WARN;
my $smooth = 100; # treePL smoothing factor
my $numsites = 670; # COI alignment size
my ( $dir, $timetree, $grafted );
GetOptions(
	'dir=s'        => \$dir, # ML_trees_collapsed
	'timetree=s'   => \$timetree, # Lepidoptera_family.nwk
	'grafted=s'    => \$grafted, # grafted.trees
	'numsites=i'   => \$numsites,
	'smooth=i'     => \$smooth,
	'verbose+'     => \$verbosity,
);

# instantiate helper objects
my $log = Bio::Phylo::Util::Logger->new(
	'-level' => $verbosity,
	'-class' => 'main',
);

$log->info("going to read time tree $timetree");
my $tt = parse_tree(
	'-format' => 'newick',
	'-file'   => $timetree,
);
$log->info("done reading time tree");

$log->info("going to read grafted tree $grafted");
my $gt = parse_tree(
	'-format' => 'newick',
	'-file'   => $grafted,
);
$log->info("done reading grafted tree");

# start iterating over families in the time tree
my @prune;
TIP: for my $tip ( @{ $tt->get_terminals } ) {
	my $family = $tip->get_name;
	my $age    = $tip->get_branch_length;
	
	# see if there is a file there
	my $file;
	if ( -e "${dir}/RAxML_bipartitions.${family}_rooted.nwk" ) {
		$file = "${dir}/RAxML_bipartitions.${family}_rooted.nwk";
		$log->info("Found file $file for family $family");
	}
	elsif ( -e "${dir}/RAxML_bipartitions.${family}.nwk" ) {
		$file = "${dir}/RAxML_bipartitions.${family}.nwk";
		$log->info("Found file $file for family $family");
	}
	else {
		push @prune, $family;
		next TIP;
	}
	
	# read the family tree file, get the ingroup, get the leftmost and rightmost tip
	my $ft = parse_tree(
		'-format' => 'newick',
		'-file'   => $file,
	);
	my ($ingroup) = sort { scalar(@{$b->get_terminals}) <=> scalar(@{$a->get_terminals}) } 
		@{ $ft->get_root->get_children };
	my $left_clade  = $ingroup->get_leftmost_terminal->get_name;
	my $right_clade = $ingroup->get_rightmost_terminal->get_name;
	$log->info("ingroup of $file is $left_clade,$right_clade");
	
	# get the equivalent nodes and MRCA in the grafted tree
	$log->info("going to look for $left_clade,$right_clade ingroup in grafted tree");	
	my $left_crown  = $gt->get_by_name( $left_clade );
	my $right_crown = $gt->get_by_name( $right_clade );
	if ( $left_crown and $right_crown ) {
		
		# annotate the MRCA
		$log->info("found MRCA ($family) in grafted tree");
		my $mrca = $gt->get_mrca([ $left_crown, $right_crown ]);
		$mrca->set_name( $family );
		$mrca->set_generic( 'age' => { $family => $age } );
	}
}

# carry over the family names
$log->info("going to carry over names");
$gt->visit_depth_first(
	'-post' => sub {
		my $node = shift;
		if ( my $fam = $node->get_generic('age') ) {
			if ( my $parent = $node->get_parent ) {
				my $pfam = $parent->get_generic('age') || {};
				for my $key ( keys %$fam ) {
					$pfam->{$key} = $fam->{$key};
				}
				$parent->set_generic( 'age' => $pfam );
			}
		}
	}
);

# find coalescent nodes that occur in the time tree
$tt->prune_tips(\@prune);
$log->info("going to find matching coalescent nodes");
my $counter = 1;
my %points;
$gt->visit_depth_first(
	'-post' => sub {
		my $node = shift;
		if ( my $fam = $node->get_generic('age') ) {
			if ( my $parent = $node->get_parent ) {
				my $pfam = $parent->get_generic('age');
				
				# we are at a node where multiple families coalesce
				if ( scalar(keys(%$fam)) < scalar(keys(%$pfam)) && ! $parent->get_generic('seen') ) {
					my @fams = keys %$pfam;
					my @tips = map { $tt->get_by_name($_) } @fams;
					my $mrca = $tt->get_mrca(\@tips);
					
					# time tree has same number of tips as coalescing families
					if ( scalar(@{$mrca->get_terminals}) == scalar(@fams) ) {
						$log->info("found equivalent clade @fams");
						my $clade = $counter++;
						my $age   = $mrca->calc_max_path_to_tips;
						my ($min) = sort { $pfam->{$a} <=> $pfam->{$b} } keys %$pfam;
						my ($max) = sort { $pfam->{$b} <=> $pfam->{$a} } keys %$pfam;
						
						$points{$clade} = {
							'left'  => $parent->get_leftmost_terminal->get_name,
							'right' => $parent->get_rightmost_terminal->get_name,
							'age'   => $age,
							'fams'  => \@fams,
							'mint'  => $min,
							'maxt'  => $max,
						};
						$parent->set_generic( 'seen' => 1 );						
					}
					else {
						$log->warn("topological conflict for @fams");
					}
				}
			}
		}		
	}
);

# scale the branch lengths
my $newick = $gt->to_newick( '-nodelabels' => 1 );

# print r8s block
print <<"R8S";
#NEXUS
BEGIN TREES;
	TREE intree = $newick
END;

BEGIN r8s;
	blformat lengths=persite nsites=${numsites} ultrametric=no;
	set smoothing=${smooth};
	collapse;
R8S

for my $p ( sort { $a <=> $b } keys %points ) {
	print "\n\t[ Coalescing families: \n";
	print "\t\t$_\n" for sort { $a cmp $b } @{ $points{$p}->{'fams'} };
	print "\t]\n";
	print "\t[ Youngest family: ", $points{$p}->{'mint'}, " ]\n";
	print "\t[ Oldest family: ", $points{$p}->{'maxt'}, " ]\n";
	print "\tMRCA clade$p ", $points{$p}->{left}, ' ', $points{$p}->{right}, ";\n";
	print "\tfixage taxon=clade$p age=", $points{$p}->{age}, ";\n";

}
print "\tdivtime method=pl algorithm=tn;\n";
print "END;\n";
