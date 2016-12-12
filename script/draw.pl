#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use Color::Spectrum;
use Bio::Phylo::Treedrawer;
use Bio::Phylo::IO qw'parse_tree';
use Bio::Phylo::Util::Logger qw':simple :levels';

# process command line arguments
my ( $intree, $fossils, $taxa );
my $verbosity = WARN;
my $width     = 2000;
my $height    = 3000;
my $padding   = 90;
my $textsize  = 18;
my $shape     = 'curvy';
my $mode      = 'phylo';
my @colors    = qw(indianred mediumpurple);
GetOptions(
	'verbose+'    => \$verbosity,
	'intree=s'    => \$intree,
	'width=i'     => \$width,
	'height=i'    => \$height,
	'padding=i'   => \$padding,
	'labelsize=i' => \$textsize,
	'shape=s'     => \$shape,
	'mode=s'      => \$mode,
	'taxa=s'      => \$taxa, # ML_trees_collapsed_families.tsv
	'fossils=s'   => \$fossils, # calibrations_timetree.tsv
	'colors=s'    => sub { @colors = split /,/, pop }
);

##########################################################################################
Bio::Phylo::Util::Logger->new(
	'-level' => $verbosity,
	'-class' => 'main',
);

INFO "going to read tree $intree";
my $tree = parse_tree(
	'-format' => 'nexus',
	'-file'   => $intree,
);

INFO "going to create tree drawer";
my $drawer = Bio::Phylo::Treedrawer->new(
	'-format'  => 'svg',
	'-width'   => $width,
	'-height'  => $height,
	'-shape'   => $shape,
	'-mode'    => $mode,
	'-tree'    => $tree,
	'-padding' => $padding,
	'-collapsed_clade_width' => 1/2,
	'-text_horiz_offset'     => 10,
	'-text_width'            => 250,
	'-branch_width'          => 4,
);

INFO "going to create color spectrum";
my $spectrum = Color::Spectrum->new;

##########################################################################################
INFO "going to set clade labels";
my $clades = 0;
my %families;
{
	open my $fh, '<', $taxa or die $!;
	while(<$fh>) {
		chomp;
		my ( $taxon, $family ) = split /\t/, $_;
		if ( my $tip = $drawer->get_tree->get_by_name( $taxon ) ) {
			if ( not $families{$family} ) {
			
				INFO "found family $family";
				$families{$family} = [];
				$clades++;
			}
			push @{ $families{$family} }, $tip;
		}
	}
	for my $fam ( keys %families ) {
		my $mrca = $drawer->get_tree->get_mrca( $families{$fam} );
		$mrca->set_clade_label($fam);
	}
}

##########################################################################################
INFO "going to paint $clades clades with spectrum (@colors)";
my @spectrum = $spectrum->generate( $clades, @colors );
my $i = 0; # incrementing index
$tree->visit_depth_first(
	'-pre' => sub {
		my $node = shift;
		if ( my $label = $node->get_clade_label ) {

			INFO "going to paint $label with $spectrum[$i]";
			$node->set_node_color($spectrum[$i]);
			for my $d ( @{ $node->get_descendants } ) {
				#$d->set_node_color($spectrum[$i]);
				$d->set_node_outline_color($spectrum[$i]);
				$d->set_branch_color($spectrum[$i]);
			}
			$node->set_font_face('Verdana');
			$node->set_font_size($textsize);
			$node->set_font_color($spectrum[$i]);
			
			
			$i++;
		}
	}
);

##########################################################################################
INFO "going to compute scale bar units";
$drawer->get_tree->ladderize('reverse');
$drawer->compute_coordinates;
my $root = $drawer->get_tree->get_root;
my ($tallest) = sort { $b->get_x <=> $a->get_x } @{ $drawer->get_tree->get_entities };
my $height_in_pixels = $tallest->get_x - $root->get_x;
my $height_in_myears = $tallest->calc_path_to_root;
my $pixels_per_myear = $height_in_pixels / $height_in_myears;

##########################################################################################
INFO "setting scale, pixels / MYA: $height_in_pixels/$height_in_myears=$pixels_per_myear";
$drawer->set_scale_options(
	'-width'   => '100%',
	'-major'   => $pixels_per_myear * 25,
	'-minor'   => $pixels_per_myear * 5,
	'-blocks'  => $pixels_per_myear * 25,
	'-label'   => 'MYA',
	'-reverse' => 1,
	'-tmpl' => sub {
	
		# round to nearest integer
		my $value = shift;
		return(($value == int($value)) ? $value : int($value + 1))
	},
	'-font' => {
		'-face'   => 'Verdana',
		'-size'   => $textsize,
		'-weight' => 'bold',
	}
);

##########################################################################################
INFO "going to collapse monophyletic genera";
$drawer->get_tree->visit_depth_first(
	'-post' => sub {
		my $node = shift;
		if ( $node->is_terminal ) {
			my $name = $node->get_name;
			if ( $name =~ /^([A-Z][a-z]+)/ ) {
				my $genus = $1;
				$node->set_generic( 'genus' => { $genus => 1 } );
				DEBUG "$name belongs to genus $genus";
			}
		}
		else {
			my %genus;
			for my $c ( @{ $node->get_children } ) {
				my $g = $c->get_generic('genus');
				for my $k ( keys %$g ) {
					$genus{$k} = $g->{$k};
				}
			}
			$node->set_generic( 'genus' => \%genus );
			DEBUG "internal node subtends: ".Dumper(\%genus);
		}
	}
);
$drawer->get_tree->visit_depth_first(
	'-pre' => sub {
		my $node = shift;
		if ( $node->is_internal ) {
			my $g = $node->get_generic('genus');
			
			# node is poly/paraphyletic
			if ( scalar(keys(%$g)) > 1 ) {
				for my $c ( @{ $node->get_children } ) {
					my $cg = $c->get_generic('genus');
					
					# child is a monophyletic internal nodes
					if ( scalar(keys(%$cg)) == 1 and $c->is_internal ) {						
						my $ntips = scalar @{ $c->get_terminals };
						my ($genus) = keys %$cg;
						my $name = "$genus (n=$ntips)";
						$c->set_name($name);
						$c->set_collapsed(1);
						$c->set_font_size($textsize);
						$c->set_font_face('Verdana');
						$c->set_font_weight('bold');
						$c->set_node_colour($c->get_node_outline_colour);
						$c->set_branch_width(0);
						INFO "found monophyletic node for $genus (n=$ntips)";
					}
				}
			}
		}
	}
);

##########################################################################################
INFO "writing SVG output to STDOUT";
for my $tip ( @{ $drawer->get_tree->get_terminals } ) {
	$tip->set_font_face('Verdana');
	$tip->set_font_size($textsize);
	$tip->set_font_style('Italic') if not $tip->get_collapsed;	
}
print $drawer->draw;
