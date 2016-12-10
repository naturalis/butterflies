#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use Color::Spectrum;
use Bio::Phylo::Treedrawer;
use Bio::Phylo::IO qw'parse_tree';
use Bio::Phylo::Util::Logger qw':simple :levels';

# process command line arguments
my ( $intree, $dir );
my $verbosity = WARN;
my $width     = 1000;
my $height    = 4000;
my $format    = 'nexus'; # of $intree
my $shape     = 'rect';
my $mode      = 'phylo';
my $prefix    = 'RAxML_bipartitions.';
my $extension = 'nwk';
my @colors    = qw(darkred darkviolet);
GetOptions(
	'verbose+'    => \$verbosity,
	'intree=s'    => \$intree,
	'dir=s'       => \$dir,
	'width=i'     => \$width,
	'height=i'    => \$height,
	'format=s'    => \$format,
	'shape=s'     => \$shape,
	'mode=s'      => \$mode,
	'prefix=s'    => \$prefix,
	'extension=s' => \$extension,
	'colors=s'    => sub { @colors = split /,/, pop }
);

# configure INFO/WARN/etc. messages 
Bio::Phylo::Util::Logger->new(
	'-level' => $verbosity,
	'-class' => 'main',
);

INFO "going to read tree $intree";
my $tree = parse_tree(
	'-format' => $format,
	'-file'   => $intree,
);

INFO "going to create tree drawer";
my $drawer = Bio::Phylo::Treedrawer->new(
	'-format' => 'svg',
	'-width'  => $width,
	'-height' => $height,
	'-shape'  => $shape,
	'-mode'   => $mode,
	'-tree'   => $tree,
);

INFO "going to create color spectrum";
my $spectrum = Color::Spectrum->new;

INFO "going to scan clades directory $dir";
my $clades = 0;
opendir my $dh, $dir or die $!;
while( my $entry = readdir $dh ) {
	if ( $entry =~ /$prefix(.+)\.$extension/ ) {
		my $family = $1;
		$family =~ s/_rooted$//g;
		if ( $family =~ /[A-Z][a-z]+/ ) {


			INFO "reading $dir/$entry for clade $family";
			my $ft = parse_tree(
				'-format' => 'newick',
				'-file'   => "$dir/$entry",
			);
			my ($ingroup) = sort {
				scalar @{ $b->get_terminals }
				<=>
				scalar @{ $a->get_terminals }
			} @{ $ft->get_root->get_children };

			INFO "matching ingroup tips with tree";
			my @matches;
			for my $tip ( @{ $ingroup->get_terminals } ) {
				my $name = $tip->get_name;
				if ( my $match = $tree->get_by_name($name) ) {

					INFO "found match for $name";
					push @matches, $match;
				}
				else {
					DEBUG "no match for $name";
				}
			}

			if ( @matches ) {

				INFO "have ingroup matches, fetching mrca";
				my $mrca = $tree->get_mrca(\@matches);
				$mrca->set_clade_label($family);
				$clades++;
			}
		}
	}
}

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
				$d->set_node_color($spectrum[$i]);
				$d->set_branch_color($spectrum[$i]);
			}
			$i++;
		}
	}
);

INFO "writing SVG output to STDOUT";
$drawer->get_tree->ladderize('reverse');
for my $tip ( @{ $drawer->get_tree->get_terminals } ) {
	$tip->set_font_face('Verdana');
	$tip->set_font_style('Italic');
	$tip->set_font_size(11);
}
print $drawer->draw;
