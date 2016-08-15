#!/usr/bin/perl
##
##
###############################################################################

use strict;

if ($#ARGV < 0) {
	print STDERR "usage: $0 <chains> <pdbfile>* \n";
	exit -1;
}

my $chnstring = shift @ARGV;
my @pdbfiles = @ARGV;

foreach my $file ( @pdbfiles ) {
	open (PDB, $file)|| die $_;
	my @buf = <PDB>;
	close (PDB);
	
	open (OUT, ">$file") || die $_;
	foreach my $line (@buf) {
		if ($line !~ /^ATOM/ && $line !~ /^HETATM/ ) {
			print OUT $line;
			next;
		}

		my $resid = substr($line, 22, 4);
		my $atmid = substr($line, 12, 4);
		my $chnid = substr($line, 21, 1);

		#if ($chnstring =~ /$chnid/) {
		if ( index($chnstring, $chnid) != -1) {
			print OUT $line;
		}
	}
}
exit 0;
