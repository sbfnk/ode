#!/usr/bin/perl

use Getopt::Long;

my $usage = "Usage: analyse_final.pl filename -a alpha -b beta -g gamma -l".
            " lambda -o omega -r rho\n";

my $filename = shift(@ARGV);
($filename eq "") && ($usage);

GetOptions("alpha=s" => \$alpha,
	   "beta=s" => \$beta,
	   "gamma=s" => \$gamma,
	   "delta=s" => \$delta,
	   "lambda=s" => \$lambda,
	   "omega=s" => \$omega,
	   "rho=s" => \$rho);

(($alpha eq "") || ($beta eq "") || ($gamma eq "") || ($delta eq "") ||
($lambda eq "") || ($omega eq "") || ($rho eq "")) && die("$usage"); 


open(IN, "<", $filename) or die "Can't open $filename for reading\n";

if ($alpha eq "v") {
  $alpha = "([0-9].[0-9][0-9])";
} else {
  $alpha = sprintf("%.2f",$alpha);
} 
if ($beta eq "v") {
  $beta = "([0-9].[0-9][0-9])";
} else {
  $beta = sprintf("%.2f",$beta);
}
if ($gamma eq "v") {
  $gamma = "([0-9].[0-9][0-9])";
} else {
  $gamma = sprintf("%.2f",$gamma);
}
if ($delta eq "v") {
  $delta = "([0-9].[0-9][0-9])";
} else {
  $delta = sprintf("%.2f",$delta);
}
if ($lambda eq "v") {
  $lambda = "([0-9].[0-9][0-9])";
} else {
  $lambda = sprintf("%.2f",$lambda);
}
if ($omega eq "v") {
  $omega = "([0-9].[0-9][0-9])";
} else {
  $omega = sprintf("%.2f",$omega);
}
if ($rho eq "v") {
  $rho = "([0-9].[0-9][0-9])";
} else {
  $rho = sprintf("%.2f",$rho);
}

my $pattern="./ab_$alpha\_$beta/dgl_$delta\_$gamma\_$lambda/or_$omega\_$rho.final\n";
my @data;

while ($dataline = <IN>) {
  $fileline = <IN>;
  if ($fileline =~ /$pattern/) {
    print "$1 $dataline";
  }
}
