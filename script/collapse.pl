#!/usr/bin/perl
use strict;
use warnings;
use File::Spec;
use Data::Dumper;
use Getopt::Long;
use LWP::UserAgent;
use List::Util 'sum';
use Bio::Phylo::Factory;
use Bio::Phylo::IO 'parse_tree';
use Bio::Phylo::Util::Logger ':levels';

# process command line arguments
my $verbosity = WARN;
my $infile;
GetOptions(
	'infile=s' => \$infile,
	'verbose+' => \$verbosity,
);

# instantiate helper objects
my $fac = Bio::Phylo::Factory->new;
my $log = Bio::Phylo::Util::Logger->new(
	'-level' => $verbosity,
	'-class' => 'main',
);

# parse input tree
$log->info("going to read newick tree from $infile");
my $tree = parse_tree(
	'-format' => 'newick',
	'-file'   => $infile,
);
$log->info("done reading newick tree from $infile");

# do monophylizer test
my %records;
$log->info("going to do remote Monophylizer test");
my $ua  = LWP::UserAgent->new;
my $res = $ua->post(
	'http://monophylizer.naturalis.nl/cgi-bin/monophylizer.pl',
	'Content_Type' => 'form-data',
	'Content' => [
		'infile' => [ 
			undef, 
			$infile,
			'Content_Type' => 'text/plain',
			'Content' => $tree->to_newick,
		],		
		'format'     => 'newick',
		'astsv'      => 'astsv',
		'submit'     => 'Go',
		'cgi'        => 1,	
		'separator'  => '|',
		'trinomials' => 1,	
	]
);
if ( $res->is_success ) {
	my @lines = split /\n/, $res->content;
	my @header;
	LINE: for my $line ( @lines ) {
		my @fields = split /\t/, $line;
		if ( not @header ) {
			@header = @fields;
			next LINE;
		}
		my %record = map { $header[$_] => $fields[$_] } 0 .. $#header;
		my $tanglees = $record{'Tanglees'};
		$record{'Tanglees'} = [ split /,/, $tanglees ];
		$records{$record{'Species'}} = \%record;
	}
	$log->info("test was successful, have results for ".scalar(keys(%records))." species");
}
else {
	$log->error("Failure: ".$res->message);
	$log->error($res->content);
	die;
}

# collect tips
$log->info("going to collect tips");
for my $tip ( @{ $tree->get_terminals } ) {
	my $label = $tip->get_name;
	$label =~ s/'//g;
	my ( $name, $id ) = split /\|/, $label;
	$name =~ s/_/ /g;
	if ( $records{$name} ) {
		$records{$name}->{'Tips'} = {} if not $records{$name}->{'Tips'};
		$records{$name}->{'Tips'}->{$id} = $tip;
	}
	else {
		$log->error("$label not seen by Monophylizer");
	}
}

# start collapsing
my %seen;
SPECIES: for my $species ( sort { scalar(@{$records{$b}->{'Tanglees'}}) <=> scalar(@{$records{$a}->{'Tanglees'}}) } keys %records ) {
	next SPECIES if $seen{$species};
	$log->info("going to collapse $species");
	
	# species is monophyletic, can just collapse as is
	my $ass = $records{$species}->{'Assessment'}; 
	if ( $ass eq 'monophyletic' ) {
		$log->info("$species is monophyletic");
		my $tips = [ values %{ $records{$species}->{'Tips'} } ];		
		
		# there is something to collapse
		if ( scalar(@$tips) > 1 ) {
			my $mrca = $tree->get_mrca([@$tips]);
			replace($species,$mrca,$tips);
		}
		
		# just rename to species
		else {
			$tips->[0]->set_name($species);
		}
		$seen{$species}++;
	}
	
	# species is para/polyphyletic, need to collect all tanglees and collapse 
	# them together
	else {
		my @tanglees = ( $species, @{ $records{$species}->{'Tanglees'} } );
		$log->info("$species is $ass, entangled with @tanglees");
		
		# collect all tips in the tangle, get their combined MRCA
		my %tips;
		for my $species ( @tanglees ) {
			$tips{$species} = [ values %{ $records{$species}->{'Tips'} } ];
		}		
		my $mrca = $tree->get_mrca([ map { @$_ } values %tips ]);
		
		# collapse each species to the MRCA
		for my $species ( keys %tips ) {
			my $tips = $tips{$species};
			replace($species,$mrca,$tips);
			$seen{$species}++;
		}
	}
}

sub replace {
	my ( $species, $mrca, $tips ) = @_;
	$log->info("$species has ".scalar(@$tips)." specimens");		

	# calculate branch length, which is the average height from the MRCA
	my @paths;
	for my $tip ( @$tips ) {
		push @paths, $tip->calc_path_to_root;
	}
	my $ptr = $mrca->calc_path_to_root;
	my $height = sum(@paths) / scalar(@paths);
	my $length = $height - $ptr;
	$log->info("average height is $height, from MRCA is $length, height of MRCA is $ptr");

	# add new child
	my $child = $fac->create_node(
		'-name' => $species,
		'-branch_length' => $length,
	);
	$mrca->set_child( $child );
	$tree->insert( $child );

	# prune old tips
	$tree->prune_tips($tips);
}

print $tree->to_newick;