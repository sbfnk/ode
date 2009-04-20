#!/usr/bin/perl

my $rho = 0.8;

my $usage = "Usage: analyse_vincent.pl filename\n";

my $filename = shift(@ARGV);
($filename ne "") || die($usage);

my $com_rho = shift(@ARGV);
($com_rho ne "") && ($rho=$com_rho+0);

open(IN, "<", $filename) or die "Can't open $filename for reading\n";

while ($line = <IN>) {
  @fields = split(' ', $line);
  my $time = 0;
  my $S = 0;
  my $I = 0;
  my $R = 0;
  my $Sprime = 0;
  my $Iprime = 0;
  my $Rprime = 0;
  my $i = 0;
  $time = shift(@fields);
  while (scalar @fields > 0) {
    my $S_i = shift(@fields);
    my $I_i = shift(@fields);
    my $R_i = shift(@fields);
    $S += $S_i;
    $I += $I_i;
    $R += $R_i;
    my $weight = 1 - ($rho**$i);
    $Sprime += $weight*$S_i;
    $Iprime += $weight*$I_i;
    $Rprime += $weight*$R_i;
    $i++;
  } 
  print "$time $S $I $R $Sprime $Iprime $Rprime\n";
}
