#!/usr/bin/perl -w

use strict;
use POSIX qw(ceil floor);

use lib (".");
use File::Basename;
use File::Path qw(rmtree);
use File::Copy qw(move copy);
use List::Util qw(max min);
use Cwd;
use Cwd 'abs_path';
use POSIX qw(ceil floor fmod fabs);
use Storable qw(dclone);
use FindBin;
use lib $FindBin::Bin;
use Getopt::Long;

my $RELAXAPP = "~/Desktop/Rosetta/main/source/bin/relax.default.linuxgccrelease".
        " -database ~/Desktop/Rosetta/main/database ";

my $WORKDIR="_temp";

#####
##### parameters
my $NMODELSOUT = 1;
my @RMS_CUTOFFS = (8,5,2.5);
my $BOLTZMANTEMP = 20;
my $REMOVEFRACTION = 0.5;

## main()
my $argc = scalar(@ARGV);
if ($argc < 1) {
    print STDERR "usage: $0 [--debug] <pdbfile>*\n";
    exit -1;
}

my ($DEBUG) = (0);
GetOptions( "debug" => \$DEBUG );

rmtree($WORKDIR);
mkdir $WORKDIR;

my @THREADED_MDLS = @ARGV;
my $nmdls = scalar (@THREADED_MDLS);
my $removeNeighbors = int($REMOVEFRACTION*$nmdls);

my %clustermap; # map each cluster to all symmetries
my (@clusterids,@symmgroups,@ddgs);
my (@modelsbyclust,@gdtsbyclust);

## stage 1 - process each cluster/symmetry seperately
my $filecounter=0;

my (%allfrags,%bbatoms,%caatoms,%resids,%fraglens,%ss);

# read in all PDBs
my $minres = 9999;
my $maxres = 0;
foreach my $pdb (@THREADED_MDLS) {
	open (PDB, $pdb) || die "Cannot open $pdb";

	$bbatoms{$pdb} = {};
	$caatoms{$pdb} = {};
	$resids{$pdb} = {};
	$allfrags{$pdb} = [];

	my $curr_frag_start = -999;
	my $last_res_read = -999;
	while (my $line = <PDB>) {
		next if ($line !~ /^ATOM/);

		my $atom = substr ($line, 12, 4);
		my $chain = substr($line, 21, 1);
		my $resid = substr($line, 17, 3);
		my $res = int( substr($line, 22, 4) );
		$minres = min($minres,$res);
		$maxres = max($maxres,$res);

		my $x = substr ($line, 30, 8);
		my $y = substr ($line, 38, 8);
		my $z = substr ($line, 46, 8);

		my $id = $atom.$chain.$res;
		$bbatoms{$pdb}->{ $id } = [$x,$y,$z];
		if ($atom eq " CA ") {
			$caatoms{$pdb}->{ $id } = [$x,$y,$z];
		}
		$resids{$pdb}->{ $id } = $resid;
		push @{$allfrags{$pdb}}, $line;
	}
	close (PDB);
}

#my $olddir = getcwd();
#chdir $WORKDIR;

# align all->all
my $npdbs = @THREADED_MDLS;
my $rmsds = []; my $nalignlen = []; my $Rs = [];
my $comis = []; my $comjs = [];
my $overlapscore = []; # for clustering

foreach my $i (0..$#THREADED_MDLS) {
foreach my $j ($i..$#THREADED_MDLS) {
	# intersection of $bbatoms{i} and {j}
	my @common_atoms = ();
	foreach ( keys %{ $caatoms{$THREADED_MDLS[$i]} } ) {
		push @common_atoms, $_ if exists $caatoms{$THREADED_MDLS[$j]}->{$_};
	}
	my $ncommon_atoms = scalar( @common_atoms );

	if ($ncommon_atoms < 4) {
		$rmsds->[$i][$j] = 9999;
		$Rs->[$i][$j] = [[1,0,0],[0,1,0],[0,0,1]];
		$nalignlen->[$i][$j] = scalar(@common_atoms);
		$overlapscore->[$i][$j] = 0;
		next;
	}

	# atom lists
	my $atoms_i = [];
	my $atoms_j = [];
	foreach ( @common_atoms ) {
		push @{ $atoms_i }, deep_copy( $caatoms{$THREADED_MDLS[$i]}->{$_} );
		push @{ $atoms_j }, deep_copy( $caatoms{$THREADED_MDLS[$j]}->{$_} );
	}

	# initial alignment
	($Rs->[$i][$j], $rmsds->[$i][$j], $comis->[$i][$j], $comjs->[$i][$j]) = rms_align( $atoms_i , $atoms_j );
	$nalignlen->[$i][$j] = scalar(@common_atoms);

	next if ($i==$j);

	# realign after trimming outliers
	foreach my $RMS_ALIGN (@RMS_CUTOFFS) {
		my @new_common_atoms = ();
		$atoms_i = [];
		$atoms_j = [];
		foreach ( @common_atoms ) {
			my $x = $caatoms{$THREADED_MDLS[$i]}->{$_};
			my $y = vadd( mapply( $Rs->[$i][$j] , vsub( $caatoms{$THREADED_MDLS[$j]}->{$_} , $comis->[$i][$j] ) ),
						  vadd($comis->[$i][$j] , $comjs->[$i][$j]) );
			my $dist = dist( $x,$y );
			if ($dist <= $RMS_ALIGN) {
				push @new_common_atoms, $_;
				push @{ $atoms_i }, deep_copy( $caatoms{$THREADED_MDLS[$i]}->{$_} );
				push @{ $atoms_j }, deep_copy( $caatoms{$THREADED_MDLS[$j]}->{$_} );
			}
		}

		$ncommon_atoms = scalar( @new_common_atoms );

		my $natoms_tot = min( scalar(keys %{ $caatoms{$THREADED_MDLS[$i]} }),  scalar(keys %{ $caatoms{$THREADED_MDLS[$j]} }) );

		if ($RMS_ALIGN == $RMS_CUTOFFS[$#RMS_CUTOFFS]) {
			$overlapscore->[$i][$j] = scalar( @new_common_atoms ) / $natoms_tot;
			last;
		}

		if ($ncommon_atoms < 4) {
			$overlapscore->[$i][$j] = 0;
			next;
		}
		@common_atoms = @new_common_atoms;

		($Rs->[$i][$j], $rmsds->[$i][$j], $comis->[$i][$j], $comjs->[$i][$j]) = rms_align( $atoms_i , $atoms_j );
		$nalignlen->[$i][$j] = scalar(@common_atoms);
	}
}
}

foreach my $i (0..$#THREADED_MDLS) {
foreach my $j (0..$i-1) {
	 $Rs->[$i][$j] =  minv( $Rs->[$j][$i] );
	 $rmsds->[$i][$j] =  deep_copy( $rmsds->[$j][$i] );
	 $comis->[$i][$j] =  vadd( $comis->[$j][$i], $comjs->[$j][$i] );
	 $comjs->[$i][$j] =  vscale( -1 , $comjs->[$j][$i] );
	 $nalignlen->[$i][$j] = $nalignlen->[$j][$i];
	 $overlapscore->[$i][$j] = $overlapscore->[$j][$i];
}
}

# start with highest probability member as the seed
my $maxprob = 1;
my $minI = 0;

## get energy-weighted median
my @towrite;
foreach my $i (0..$#THREADED_MDLS) { $towrite[$i] = 0; }

my @chooseModel = (-1) x $NMODELSOUT;
my $modelswritten = 0;

while ($modelswritten < $NMODELSOUT) {
	last if ($modelswritten == $NMODELSOUT);

	my $lowmaxsub = -999999;
	outer: foreach my $i (0..$#THREADED_MDLS) {
		next if ($towrite[$i] != 0); # already accounted for

		my $pdbidstem = $THREADED_MDLS[$i]; $pdbidstem =~ s/\.pdb$//g;
		my $score = 1;
		foreach my $j (0..$#THREADED_MDLS) {
			next if ($towrite[$j] != 0);
			next if ($i==$j); # already counted this

			my $pdbidstemJ = $THREADED_MDLS[$j];
			$pdbidstemJ =~ s/\.pdb$//g;
			$score += $overlapscore->[$i][$j];
		}

		if ($score > $lowmaxsub) {
			$chooseModel[$modelswritten] = $i;
			$lowmaxsub = $score;
		}
	}
	$towrite[$chooseModel[$modelswritten]] = $modelswritten+1;

	# remove this model and its neighbors
	my $tag = $THREADED_MDLS[$chooseModel[$modelswritten]];
	$tag =~ s/\.pdb$//g;
	my %to_remove = ($tag=>1);

	my @all_outfiles = ("$tag.pdb");
	removeloop: foreach my $k (1..$removeNeighbors) {
		my ($maxScore,$maxJ) = (0,-1);
		foreach my $j (0..$#THREADED_MDLS) {
			next if ($towrite[$j] != 0);
			next if ($chooseModel[$modelswritten]==$j);
			my $score = $overlapscore->[$chooseModel[$modelswritten]][$j];
			if ($score > $maxScore) { $maxScore = $score; $maxJ = $j; }
		}
		if ($maxJ != -1) {
			$towrite[$maxJ] = -1;
			$tag = $THREADED_MDLS[$maxJ];

			# align and write
			# add j to alignment
			my $R = $Rs->[$chooseModel[$modelswritten]][$maxJ];
			my $preT = $comis->[$chooseModel[$modelswritten]][$maxJ];
			my $postT = vadd( $comjs->[$chooseModel[$modelswritten]][$maxJ], $preT );

			open INPDB, "$tag";
			my @pdblines = <INPDB>;
			close INPDB;

			my $outfile = "$WORKDIR/model".($modelswritten+1)."_".$k.".pdb";
			push @all_outfiles , $outfile;
			open OUTPDB, ">$outfile";
			foreach my $line (@pdblines) {
				if ($line =~ /^ATOM/) {
					my $newX = [ substr ($line, 30, 8), substr ($line, 38, 8), substr ($line, 46, 8) ];
					$newX = vsub( $newX, $preT);
					$newX = mapply( $R, $newX );
					$newX = vadd( $newX, $postT);
					substr ($line, 30, 8) = sprintf ("%8.3f", $newX->[0]);
					substr ($line, 38, 8) = sprintf ("%8.3f", $newX->[1]);
					substr ($line, 46, 8) = sprintf ("%8.3f", $newX->[2]);
				}
				print OUTPDB $line."\n";
			}
			close OUTPDB;

			$tag =~ s/\.pdb$//g;
			$to_remove{$tag} = 1;
		}
	}

	my $bfact_mdl = "$WORKDIR/bfact".($modelswritten+1).".pdb";
	my $avg_mdl = "average".($modelswritten+1).".pdb";
	average_model($avg_mdl,$bfact_mdl , \@all_outfiles);

	$modelswritten++;
}

$modelsbyclust[$filecounter] = $WORKDIR;
$filecounter++;
#chdir $olddir;

# relax
cart_relax("average1.pdb");

if (!$DEBUG) { rmtree($WORKDIR); }

exit 0;

###########
###########
###########
###########
###########
###########

sub dist {
    my ($x, $y) = @_;
    my $z = [ $x->[0]-$y->[0] , $x->[1]-$y->[1] , $x->[2]-$y->[2] ];
    return sqrt( $z->[0]*$z->[0] + $z->[1]*$z->[1] + $z->[2]*$z->[2] );
}

sub cart_relax {
    my ($infile) = @_;
    my $cmd = "$RELAXAPP -in:file:s $infile -relax:cartesian -score:weights talaris2013 -set_weights cart_bonded 0.5 pro_close 0 -default_max_cycles 200 -default_repeats 1 -crystal_refine -relax::constrain_relax_to_start_coords -relax::ramp_constraints false";
    `$cmd`;
}


sub average_model {
	my ($outmodel, $outmodelB, $inmodels) = @_;
	my %models;
	my %lines;
	my @idlist;
	my ($minres,$maxres) = (9999,0);

	my $firstmodel=1;
	foreach my $pdbfile (@{ $inmodels }) {
		open (PDB, $pdbfile) || print STDERR "Cannot open $pdbfile";
		while (<PDB>) {
			if (! /^ATOM/) {
				next;
			}
			my $atom = substr ($_, 13, 3);
			my $chain = substr($_, 21, 1);
			my $res = int( substr($_, 22, 4) );
			my $id = $res.$atom;
			$minres = min( $minres, $res );
			$maxres = max( $maxres, $res );

			my $x = substr ($_, 30, 8);
			my $y = substr ($_, 38, 8);
			my $z = substr ($_, 46, 8);

			push @{ $models{ $id } }, [$x,$y,$z];
			if ($firstmodel == 1) { $lines{ $id } = $_; };
		}
		$firstmodel=0;
		close (PDB);
	}

	my %Bs;
	open (OUT1, ">$outmodel") || print STDERR "Cannot open $_";
	foreach my $resid ( $minres .. $maxres ) {
		foreach my $atomid ("N  ","CA ","C  ","O  ") {
			next if (!defined $models{ $resid.$atomid });

			my $sumX=[0,0,0];
			my $sumX2=[0,0,0];
			my $N = scalar @{ $models{ $resid.$atomid } };
			next if ($N<1);

			foreach my $coord (@{ $models{ $resid.$atomid } }) {
				foreach my $i (0..2) {
					$sumX->[$i] += $coord->[$i];
					$sumX2->[$i] += $coord->[$i]*$coord->[$i];
				}
			}
			foreach my $i (0..2) {
				$sumX->[$i] /= $N;
				$sumX2->[$i] = $sumX2->[$i]/$N - $sumX->[$i]*$sumX->[$i];
			}
			my $B = 8/3 * (3.1415) * (3.1415) * ($sumX2->[0] + $sumX2->[1] + $sumX2->[2]);

			# if there is too much diversity, just use model 1
			if ($B>50.0) {
				$sumX = $models{ $resid.$atomid }->[0];
			}
			if ($B>400.0) { $B=400; }
			$Bs{$resid.$atomid} = $B;

			my $average_line = $lines{ $resid.$atomid };
			my $Bfact_line = $lines{ $resid.$atomid };

			substr ($average_line, 30, 8) = sprintf ("%8.3f", $sumX->[0]);
			substr ($average_line, 38, 8) = sprintf ("%8.3f", $sumX->[1]);
			substr ($average_line, 46, 8) = sprintf ("%8.3f", $sumX->[2]);
			substr ($average_line, 60, 7) = sprintf ("%7.3f", $B);

			print OUT1 $average_line."\n";
		}
	}
	close (OUT1);

	open (OUT2, ">$outmodelB") || print STDERR "Cannot open $_";
	open (PDB, $inmodels->[0]) || print STDERR "Cannot open ".$inmodels->[0];
		while (<PDB>) {
		if (! /^ATOM/) {
			next;
		}

		my $atom = substr ($_, 13, 3);
		my $eff_atom = $atom;
		my $multiplier = 1.0;

		if ($atom ne "N  " && $atom ne "CA " && $atom ne "C  " && $atom ne "O  ") {
			$eff_atom = "CA ";
			if (substr($atom,1,1) eq 'B') {$multiplier = 1.1;}
			elsif (substr($atom,1,1) eq 'G') {$multiplier = 1.2;}
			else {$multiplier = 2.0;}
			if (substr($atom,0,1) eq 'H') {$multiplier *= 1.2;}
		}

		my $res = int( substr($_, 22, 4) );

		my $B = $Bs{$res.$eff_atom}*$multiplier;
		if ($B>400.0) { $B=400; }

		substr ($_, 60, 7) = sprintf ("%7.3f", $B);
		print OUT2 $_;
	}

	close (OUT2);
}



sub deep_copy {
	my $this = shift;
	if (not ref $this) {
		$this;
	} elsif (ref $this eq "ARRAY") {
		[map deep_copy($_), @$this];
	} elsif (ref $this eq "HASH") {
		+{map { $_ => deep_copy($this->{$_}) } keys %$this};
	} else { die "what type is $_?" }
}


# rotation from euler angles
sub euler {
	my ($aa, $bb, $gg) = @_;
	my $MM;

	$MM->[0][0] = (-sin($aa)*cos($bb)*sin($gg) + cos($aa)*cos($gg));
	$MM->[0][1] = ( cos($aa)*cos($bb)*sin($gg) + sin($aa)*cos($gg));
	$MM->[0][2] = ( sin($bb)*sin($gg));
	$MM->[1][0] = (-sin($aa)*cos($bb)*cos($gg) - cos($aa)*sin($gg));
	$MM->[1][1] = ( cos($aa)*cos($bb)*cos($gg) - sin($aa)*sin($gg));
	$MM->[1][2] = ( sin($bb)*cos($gg));
	$MM->[2][0] = ( sin($aa)*sin($bb));
	$MM->[2][1] = (-cos($aa)*sin($bb));
	$MM->[2][2] = ( cos($bb));

	return $MM;
}

# my ($X,$Y,$Z,$W)=R2quat($M)
sub R2quat {
	my $R = shift;
	my ($S,$X,$Y,$Z,$W);
	if ( $R->[0][0] > $R->[1][1] && $R->[0][0] > $R->[2][2] )  {
		$S  = sqrt( 1.0 + $R->[0][0] - $R->[1][1] - $R->[2][2] ) * 2;
		$X = 0.25 * $S;
		$Y = ($R->[1][0] + $R->[0][1] ) / $S;
		$Z = ($R->[2][0] + $R->[0][2] ) / $S;
		$W = ($R->[2][1] - $R->[1][2] ) / $S;
	} elsif ( $R->[1][1] > $R->[2][2] ) {
		$S  = sqrt( 1.0 + $R->[1][1] - $R->[0][0] - $R->[2][2] ) * 2;
		$X = ($R->[1][0] + $R->[0][1] ) / $S;
		$Y = 0.25 * $S;
		$Z = ($R->[2][1] + $R->[1][2] ) / $S;
		$W = ($R->[0][2] - $R->[2][0] ) / $S;
	} else {
		$S  = sqrt( 1.0 + $R->[2][2] - $R->[0][0] - $R->[1][1] ) * 2;
		$X = ($R->[0][2] + $R->[2][0] ) / $S;
		$Y = ($R->[2][1] + $R->[1][2] ) / $S;
		$Z = 0.25 * $S;
		$W = ($R->[1][0] - $R->[0][1]) / $S;
	}
	return ($X,$Y,$Z,$W);
}


sub R2angle {
	my $R = shift;
	my $x = $R->[2][1]-$R->[1][2];
	my $y = $R->[0][2]-$R->[2][0];
	my $z = $R->[1][0]-$R->[0][1];
	my $r = sqrt($x*$x + $y*$y + $z*$z);
	my $t = $R->[0][0]+$R->[1][1]+$R->[2][2];
	my $th=atan2($r,$t-1);
	return $th;
}


# my ($R)=R2quat($X,$Y,$Z,$W)
sub quat2R {
	my ($X,$Y,$Z,$W) = @_;
	my $xx = $X * $X; my $xy = $X * $Y; my $xz = $X * $Z;
	my $xw = $X * $W; my $yy = $Y * $Y; my $yz = $Y * $Z;
	my $yw = $Y * $W; my $zz = $Z * $Z; my $zw = $Z * $W;
	my $R = [ [ 1 - 2 * ( $yy+$zz ) ,     2 * ( $xy-$zw ) ,     2 * ( $xz+$yw ) ] ,
	          [     2 * ( $xy+$zw ) , 1 - 2 * ( $xx+$zz ) ,     2 * ( $yz-$xw ) ] ,
	          [     2 * ( $xz-$yw ) ,     2 * ( $yz+$xw ) , 1 - 2 * ( $xx+$yy ) ] ];
	return $R;
}

# my ($R)=R2quat($X,$Y,$Z,$W)
sub quatnorm {
	my ($X,$Y,$Z,$W) = @_;
	my $S = sqrt( $X*$X+$Y*$Y+$Z*$Z+$W*$W );
	return [ $X/$S , $Y/$S , $Z/$S , $W/$S ];
}

#####################################
#####################################

# vector addition
sub vadd {
	my ($x, $y) = @_;
	return [ $x->[0]+$y->[0], $x->[1]+$y->[1], $x->[2]+$y->[2] ];
}

# vector subtraction
sub vsub {
	my ($x, $y) = @_;
	return [ $x->[0]-$y->[0], $x->[1]-$y->[1], $x->[2]-$y->[2] ];
}

# mult vector by scalar
sub vscale {
	my ($x, $y) = @_;
	return [ $x*$y->[0], $x*$y->[1], $x*$y->[2] ];
}

# "min mod"
sub minmod {
	my ($x,$y) = @_;
	my $r = fmod($x,$y);
	if ($r < -fabs( $y/2.0 ) ) { $r += fabs( $y ); }
	elsif ($r >  fabs( $y/2.0 ) ) { $r -= fabs( $y ); }
	return $r;
}

# vector min-modulus
sub vminmod {
	my ($x,$y) = @_;
	return [ minmod($x->[0],$y->[0]), minmod($x->[1],$y->[1]), minmod($x->[2],$y->[2]) ];
}


#####################################
#####################################

# raise a matrix to a power
# dumb way of doing it
sub mpow {
	my ($mat, $pow) = @_;
	my $matpow = $mat;
	foreach my $i (2..$pow) {
		$matpow = mmult( $mat, $matpow );
	}
	return $matpow;
}

# matrix x vector mult
sub mapply {
	my ($rotmat, $cart) = @_;
	my $out = [0, 0, 0];
	my ($i, $j);
	for ($i=0; $i < 3; ++$i) {
		for ($j=0; $j < 3; ++$j) {
			$out->[$i] += $rotmat->[$i][$j] * $cart->[$j];
		}
	}
	return $out;
}

# matrix x matrix mult
sub mmult {
	my ($m1, $m2) = @_;
	my $out = [ [0,0,0], [0,0,0], [0,0,0] ];
	my ($i, $j, $k);
	for ($i=0; $i<3; ++$i) {
		for ($j=0; $j<3; ++$j) {
			for ($k=0; $k<3; ++$k) {
				$out->[$i][$j] += $m1->[$i][$k] * $m2->[$k][$j];
			}
		}
	}
	return $out;
}



# matrix inversion
sub minv {
	my $M = shift;
	my $Minv = [ [1,0,0] , [0,1,0] , [0,0,1] ];
	my $D = $M->[0][0] * ( $M->[1][1]*$M->[2][2] - $M->[2][1]*$M->[1][2] ) -
		    $M->[0][1] * ( $M->[1][0]*$M->[2][2] - $M->[1][2]*$M->[2][0] ) +
		    $M->[0][2] * ( $M->[1][0]*$M->[2][1] - $M->[1][1]*$M->[2][0] );
	if ($D == 0)  {
		print STDERR "ERROR ... Inversion of singular matrix!\n";
		exit -1;
	}

	$Minv->[0][0] =  ($M->[1][1]*$M->[2][2]-$M->[1][2]*$M->[2][1])/$D;
	$Minv->[0][1] = -($M->[0][1]*$M->[2][2]-$M->[0][2]*$M->[2][1])/$D;
	$Minv->[0][2] =  ($M->[0][1]*$M->[1][2]-$M->[0][2]*$M->[1][1])/$D;
	$Minv->[1][0] = -($M->[1][0]*$M->[2][2]-$M->[2][0]*$M->[1][2])/$D;
	$Minv->[1][1] =  ($M->[0][0]*$M->[2][2]-$M->[0][2]*$M->[2][0])/$D;
	$Minv->[1][2] = -($M->[0][0]*$M->[1][2]-$M->[0][2]*$M->[1][0])/$D;
	$Minv->[2][0] =  ($M->[1][0]*$M->[2][1]-$M->[2][0]*$M->[1][1])/$D;
	$Minv->[2][1] = -($M->[0][0]*$M->[2][1]-$M->[0][1]*$M->[2][0])/$D;
	$Minv->[2][2] =  ($M->[0][0]*$M->[1][1]-$M->[0][1]*$M->[1][0])/$D;

	return $Minv;
}

# vector norm
sub vnorm {
	my $x = shift;
	return sqrt( ($x->[0]*$x->[0]) + ($x->[1]*$x->[1]) + ($x->[2]*$x->[2]) );
}

# vector norm^2
sub vnorm2 {
	my $x = shift;
	return ( ($x->[0]*$x->[0]) + ($x->[1]*$x->[1]) + ($x->[2]*$x->[2]) );
}

# cart distance
sub vdist {
	my ($x1, $x2) = @_;
	return sqrt (
	           ($x1->[0]-$x2->[0])*($x1->[0]-$x2->[0]) +
	           ($x1->[1]-$x2->[1])*($x1->[1]-$x2->[1]) +
	           ($x1->[2]-$x2->[2])*($x1->[2]-$x2->[2])
	            );
}

# are two transformations the inverso of one another
sub is_inverse {
	my $tol = 1e-8;
	my ( $R_i,$T_i, $R_j,$T_j ) = @_;

	my $testR = mmult( $R_i , $R_j );
	my $testT = vadd( mapply( $R_j,$T_i ) , $T_j );

	my $errR = square($testR->[0][0]-1) + square($testR->[1][1]-1) + square($testR->[2][2]-1) +
	           square($testR->[0][1])   + square($testR->[0][2])   + square($testR->[1][2]) +
	           square($testR->[1][0])   + square($testR->[2][0])   + square($testR->[2][1]);
	my $errT = square($testT->[0])+square($testT->[1])+square($testT->[2]);

#print " (0) $errR   $errT\n";
	if ($errR < $tol && $errT < $tol) { return 1; }
	return 0;
}

#
sub is_identity {
	my $testR = shift;
	my $tol = 1e-3;
	if (scalar( @_ ) >= 1) {
		$tol = shift;
	}
	my $errR = square($testR->[0][0]-1) + square($testR->[1][1]-1) + square($testR->[2][2]-1) +
	           square($testR->[0][1])   + square($testR->[0][2])   + square($testR->[1][2]) +
	           square($testR->[1][0])   + square($testR->[2][0])   + square($testR->[2][1]);

	if ($errR < $tol) { return 1; }
	return 0;
}

# is the transform (Rn,Tn) equivalent to the transform (Ri,Ti)->(Rj,Tj)
sub is_equivalent {
	my $tol = 1e-8;
	my ( $R_n,$T_n, $R_i,$T_i, $R_j,$T_j ) = @_;

	my $R_i_inv = minv( $R_i );
	my $T_i_inv = [ -$T_i->[0], -$T_i->[1], -$T_i->[2] ];

	my $R_test = mmult( $R_i_inv, $R_j );
	my $T_test = mapply( $R_i_inv, vsub( $T_j, $T_i ) );

	my $errR = square($R_test->[0][0]-$R_n->[0][0])+square($R_test->[0][1]-$R_n->[0][1])+square($R_test->[0][2]-$R_n->[0][2])
	         + square($R_test->[1][0]-$R_n->[1][0])+square($R_test->[1][1]-$R_n->[1][1])+square($R_test->[1][2]-$R_n->[1][2])
	         + square($R_test->[2][0]-$R_n->[2][0])+square($R_test->[2][1]-$R_n->[2][1])+square($R_test->[2][2]-$R_n->[2][2]);
	my $errT = square($T_test->[0]-$T_n->[0])+square($T_test->[1]-$T_n->[1])+square($T_test->[2]-$T_n->[2]);

	##
	## don't think we need to check T ...
	if ($errR < $tol) { return 1; }
	#if ($errR < $tol && $errT < $tol) { return 1; }
	return 0;
}



# my ($R,$rmsd, $Ycom, $Ycom_to_Xcom) = rms_align( $x,$y );
sub rms_align {
	my ($X,$Y) = @_;

	my ($nlist,$mov_com, $mov_to_ref, $R, $E0) = setup_rotation( $X, $Y );
	my ($U,$residual) = calculate_rotation_matrix($R,$E0);

	my $rmsd = sqrt( fabs($residual*2.0/$nlist) );
	return ($U,$rmsd,$mov_com, $mov_to_ref);
}

# my $com = recenter( $x );
sub recenter {
	my ($X) = @_;
	my $Natms = scalar(@{ $X });
	my $com = [0,0,0];

	foreach my $n (0..$Natms-1) {
		foreach my $i (0..2) {
			$com->[$i] += $X->[$n][$i];
		}
	}

	foreach my $i (0..2) {
		$com->[$i] /= $Natms;
	}

	foreach my $n (0..$Natms-1) {
		foreach my $i (0..2) {
			$X->[$n][$i] -= $com->[$i];
		}
	}
	return $com;
}


# normalize($a)
sub normalize {
	my $a = shift;
	my $b = sqrt($a->[0]*$a->[0] + $a->[1]*$a->[1] + $a->[2]*$a->[2]);
	if ($b > 1e-6) {
		$a->[0] /= $b; $a->[1] /= $b; $a->[2] /= $b;
	}
}


# my $a_dot_b = dot($a,$b)
sub dot {
	my ($a,$b) = @_;
	return ($a->[0]*$b->[0] + $a->[1]*$b->[1] + $a->[2]*$b->[2]);
}


# my $a = cross ( b , c )
sub cross {
	my ($b,$c) = @_;
	my $a = [ $b->[1]*$c->[2] - $b->[2]*$c->[1] ,
	          $b->[2]*$c->[0] - $b->[0]*$c->[2] ,
	          $b->[0]*$c->[1] - $b->[1]*$c->[0] ];
	return $a;
}



# ($nlist,$mov_com, $mov_to_ref, $R, $E0) = setup_rotation( $ref_xlist, $mov_xlist )
sub setup_rotation {
	my ( $ref_xlist, $mov_xlist ) = @_;

	my $nlist = min( scalar(@{ $ref_xlist }) , scalar(@{ $mov_xlist }) );
	my $ref_com = [0,0,0];
	my $mov_com = [0,0,0];
	my $mov_to_ref = [0,0,0];

	foreach my $n (0..$nlist-1) {
		foreach my $i (0..2) {
			$mov_com->[$i] += $mov_xlist->[$n][$i];
			$ref_com->[$i] += $ref_xlist->[$n][$i];
		}
    }
	foreach my $i (0..2) {
		$mov_com->[$i] /= $nlist;
		$ref_com->[$i] /= $nlist;
		$mov_to_ref->[$i] = $ref_com->[$i] - $mov_com->[$i];
	}


	# shift mov_xlist and ref_xlist to centre of mass */
	foreach my $n (0..$nlist-1) {
		foreach my $i (0..2) {
			$mov_xlist->[$n][$i] -= $mov_com->[$i];
			$ref_xlist->[$n][$i] -= $ref_com->[$i];
		}
    }

	# initialize
	my $R = [ [0,0,0] , [0,0,0] , [0,0,0] ];
	my $E0 = 0.0;

	foreach my $n (0..$nlist-1) {
		# E0 = 1/2 * sum(over n): y(n)*y(n) + x(n)*x(n)
		foreach my $i (0..2) {
		  $E0 +=  $mov_xlist->[$n][$i] * $mov_xlist->[$n][$i]
				+ $ref_xlist->[$n][$i] * $ref_xlist->[$n][$i];
		}

		# R[i,j] = sum(over n): y(n,i) * x(n,j)
		foreach my $i (0..2) {
			foreach my $j (0..2) {
			   $R->[$i][$j] += $mov_xlist->[$n][$i] * $ref_xlist->[$n][$j];
			}
		}
	}
	$E0 *= 0.5;

	return ($nlist,$mov_com, $mov_to_ref, $R, $E0);
}

# helper funct
sub j_rotate {
	my ($a,$i,$j,$k,$l,$s,$tau) = @_;
	my $g = $a->[$i][$j];
	my $h = $a->[$k][$l];
	$a->[$i][$j] = $g-$s*($h+$g*$tau);
    $a->[$k][$l] = $h+$s*($g-$h*$tau);
}

# ($d,$v,$nrot) = jacobi3($a)
#    computes eigenval and eigen_vec of a real 3x3
#    symmetric matrix. On output, elements of a that are above
#    the diagonal are destroyed. d[1..3] returns the
#    eigenval of a. v[1..3][1..3] is a matrix whose
#    columns contain, on output, the normalized eigen_vec of a.
#    n_rot returns the number of Jacobi rotations that were required
sub jacobi3 {
	my $a = shift;

	my $v = [ [1,0,0] , [0,1,0] , [0,0,1] ];
	my $b = [ $a->[0][0] , $a->[1][1] , $a->[2][2] ];
	my $d = [ $a->[0][0] , $a->[1][1] , $a->[2][2] ];
	my $z = [0,0,0];
	my $n_rot = 0;
	my $thresh = 0;

	# 50 tries!
	foreach my $count (0..49) {

		# sum off-diagonal elements
		my $sum = fabs($a->[0][1])+fabs($a->[0][2])+fabs($a->[1][2]);

		# if converged to machine underflow
		if ($sum == 0.0) {
	      return($d,$v,$n_rot);
		}

		# on 1st three sweeps..
		my $thresh = 0;
		if ($count < 3) {
			$thresh = $sum * 0.2 / 9.0;
		}

		foreach my $i (0,1) {
			foreach my $j ($i+1..2) {
				my $g = 100.0 * fabs($a->[$i][$j]);

				# after four sweeps, skip the rotation if
				# the off-diagonal element is small
				if ( $count > 3
				      && fabs($d->[$i])+$g == fabs($d->[$i])
				      && fabs($d->[$j])+$g == fabs($d->[$j]) ) {
					$a->[$i][$j] = 0.0;
				} elsif (fabs($a->[$i][$j]) > $thresh) {
					my $h = $d->[$j] - $d->[$i];
					my ($t,$s,$tau,$theta);

					if (fabs($h)+$g == fabs($h)) {
						$t = $a->[$i][$j] / $h;
					} else {
						$theta = 0.5 * $h / ($a->[$i][$j]);
						$t = 1.0 / ( fabs($theta) + sqrt(1.0 + $theta*$theta) );
						if ($theta < 0.0) { $t = -$t; }
					}

					my $c = 1.0 / sqrt(1 + $t*$t);
					$s = $t * $c;
					$tau = $s / (1.0 + $c);
					$h = $t * $a->[$i][$j];

					$z->[$i] -= $h;
					$z->[$j] += $h;
					$d->[$i] -= $h;
					$d->[$j] += $h;

					$a->[$i][$j] = 0.0;

					foreach my $k (0..$i-1) {
						j_rotate($a, $k, $i, $k, $j, $s, $tau);
					}
					foreach my $k ($i+1..$j-1) {
						j_rotate($a, $i, $k, $k, $j, $s, $tau);
					}
					foreach my $k ($j+1..2) {
						j_rotate($a, $i, $k, $j, $k, $s, $tau);
					}
					foreach my $k (0..2) {
						j_rotate($v, $k, $i, $k, $j, $s, $tau);
					}
					$n_rot++;
				}
			}
		}

		foreach my $i (0..2) {
			$b->[$i] += $z->[$i];
			$d->[$i] = $b->[$i];
			$z->[$i] = 0.0;
		}
	}

	print STDERR "WARNING: Too many iterations in jacobi3!  You're bad and you should feel bad.\n";
	exit -1;
}



# ($eigen_vec, $eigenval) = diagonalize_symmetric( $matrix )
sub diagonalize_symmetric {
	my $matrix = shift;
	my ($eigenval,$vec,$n_rot) = jacobi3($matrix);

	# sort solutions by eigenval
	foreach my $i (0..2) {
		my $k = $i;
		my $val = $eigenval->[$i];

		foreach my $j ($i+1..2) {
			if ($eigenval->[$j] >= $val) {
				$k = $j;
				$val = $eigenval->[$k];
			}
		}

		if ($k != $i) {
			$eigenval->[$k] = $eigenval->[$i];
			$eigenval->[$i] = $val;
			foreach my $j (0..2) {
				$val = $vec->[$j][$i];
				$vec->[$j][$i] = $vec->[$j][$k];
				$vec->[$j][$k] = $val;
			}
		}
	}

	# transpose
	my $eigen_vec = [ [$vec->[0][0],$vec->[1][0],$vec->[2][0]] ,
	                  [$vec->[0][1],$vec->[1][1],$vec->[2][1]] ,
	                  [$vec->[0][2],$vec->[1][2],$vec->[2][2]] ];
	return ($eigen_vec, $eigenval);
}


# ($U,$residual) = calculate_rotation_matrix($R,$E0)
sub calculate_rotation_matrix {
	my ($R,$E0) = @_;

	my $RtR = [ [0,0,0] , [0,0,0] , [0,0,0] ];
	my $left_eigenvec = [ [0,0,0] , [0,0,0] , [0,0,0] ];
	my $right_eigenvec = [ [0,0,0] , [0,0,0] , [0,0,0] ];
	my $eigenval = 0;

	 # Rt <- transpose of R
	my $Rt = [ [$R->[0][0],$R->[1][0],$R->[2][0]] ,
	           [$R->[0][1],$R->[1][1],$R->[2][1]] ,
	           [$R->[0][2],$R->[1][2],$R->[2][2]] ];

	# make symmetric RtR = Rt X R
	foreach my $i (0..2) {
		foreach my $j (0..2) {
			$RtR->[$i][$j] = 0.0;
			foreach my $k (0..2) {
				$RtR->[$i][$j] += $Rt->[$k][$i] * $R->[$j][$k];
			}
		}
	}

	($right_eigenvec, $eigenval) = diagonalize_symmetric( $RtR );

	# right_eigenvec's should be an orthogonal system but could be left
	#   or right-handed. Let's force into right-handed system.
	$right_eigenvec->[2] = cross($right_eigenvec->[0], $right_eigenvec->[1]);

	# From the Kabsch algorithm, the eigenvec's of RtR
	#   are identical to the right_eigenvec's of R.
	#   This means that left_eigenvec = R x right_eigenvec
	foreach my $i (0..2) {
		foreach my $j (0..2) {
			$left_eigenvec->[$i][$j] = dot($right_eigenvec->[$i], $Rt->[$j]);
		}
	}

	foreach my $i (0..2) {
		normalize($left_eigenvec->[$i]);
	}

	# Force left_eigenvec[2] to be orthogonal to the other vectors.
	# First check if the rotational matrices generated from the
	#   orthogonal eigenvectors are in a right-handed or left-handed
	#   coordinate system - given by sigma. Sigma is needed to
	#   resolve this ambiguity in calculating the RMSD.
	my $sigma = 1.0;
	my $v = cross($left_eigenvec->[0], $left_eigenvec->[1]);
	if (dot($v, $left_eigenvec->[2]) < 0.0) {
		$sigma = -1.0;
	}
	foreach my $i (0..2) {
	    $left_eigenvec->[2][$i] = $v->[$i];
	}

	# calc optimal rotation matrix U that minimises residual
	my $U = [ [0,0,0] , [0,0,0] , [0,0,0] ];
	foreach my $i (0..2) {
		foreach my $j (0..2) {
			foreach my $k (0..2) {
				$U->[$i][$j] += $left_eigenvec->[$k][$i] * $right_eigenvec->[$k][$j];
    		}
		}
	}

	my $residual = $E0 - sqrt(fabs($eigenval->[0]))
	                   - sqrt(fabs($eigenval->[1]))
	                   - $sigma * sqrt(fabs($eigenval->[2]));

	return ($U,$residual);
}

