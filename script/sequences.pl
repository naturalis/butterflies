#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;

# process command line arguments
my ( $table, $sequences );
GetOptions(
	'table=s'     => \$table,
	'sequences=s' => \$sequences,
);

# read table
my %exemplars;
{
	open my $fh, '<', $table or die $!;
	while(<$fh>) {
		chomp;
		my ( $tree, $taxon ) = split /\t/, $_;
		$taxon =~ s/_/ /g;
		$exemplars{$taxon} = { 'tree' => $tree };
	}
}

# read fasta
open my $fh, '<', $sequences or die $!;
my ( $defline, $seq );
while(<$fh>) {
	if ( /^(>.+)$/ ) {
		my $newline = $1;
		if ( $defline and $seq ) {
			my @parts = split /\|/, $defline;
			if ( $exemplars{$parts[1]} ) {
				perhaps_store( $exemplars{$parts[1]}, $seq, $defline );
			}
		}
		$defline = $newline;
		$seq = '';
	}
	elsif ( eof $fh ) {
		if ( $defline and $seq ) {
			my @parts = split /\|/, $defline;
			if ( $exemplars{$parts[1]} ) {
				perhaps_store( $exemplars{$parts[1]}, $seq, $defline );
			}
		}	
	}
	else {
		$seq .= $_;
	}
}

sub perhaps_store {
	my ( $obj, $seq, $defline ) = @_;
	if ( not $obj->{'seq'} ) {
		$obj->{'seq'} = $seq;
		$obj->{'def'} = $defline;
	}
	else {
		my $current_gaps = $obj->{'seq'} =~ tr/-/-/;
		$current_gaps += $obj->{'seq'} =~ tr/N/N/;		
		my $future_gaps = $seq =~ tr/-/-/;
		$future_gaps += $seq =~ tr/N/N/;
		if ( $future_gaps < $current_gaps ) {
			$obj->{'seq'} = $seq;
			$obj->{'def'} = $defline;		
		}
	}
}

for my $obj ( values %exemplars ) {
	if ( not $obj->{'seq'} ) {
		die;
	}
	print $obj->{'def'}, "\n", $obj->{'seq'};
}