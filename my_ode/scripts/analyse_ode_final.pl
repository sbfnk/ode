#!/usr/bin/perl

my $usage = "Usage: analyse_vincent.pl filename\n";

my $filename = shift(@ARGV);
($filename ne "") || die($usage);

my $com_rho = shift(@ARGV);
if ($com_rho ne "") {
  $rho=$com_rho+0;
} else {
  $rho = 0;
}

open(IN, "<", $filename) or die "Can't open $filename for reading\n";

while ($line = <IN>) {
  @fields = split(' ', $line);
  my $S = 0;
  my $I = 0;
  my $R = 0;
  my $Sprime = 0;
  my $Iprime = 0;
  my $Rprime = 0;
  my $i = 0;
  while (scalar @fields > 0) {
    my $S_i = shift(@fields);
    my $I_i = shift(@fields);
    my $R_i = shift(@fields);
    $S += $S_i;
    $I += $I_i;
    $R += $R_i;
    if ($rho > 0) {
      my $weight = 1 - ($rho**$i);
      $Sprime += $weight*$S_i;
      $Iprime += $weight*$I_i;
      $Rprime += $weight*$R_i;
      $i++;
    } 
  }
  print "$S $I $R";
  if ($rho > 0) {
    print " $Sprime $Iprime $Rprime";
  }
  print "\n";
}
