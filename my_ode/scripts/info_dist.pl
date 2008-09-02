#!/usr/bin/perl

my $usage = "Usage: analyse_vincent.pl filename time\n";

my $filename = shift(@ARGV);
($filename ne "") || die($usage);

my $time = (shift(@ARGV))+0;

open(IN, "<", $filename) or die "Can't open $filename for reading\n";

my @fields;
my $curtime = $time - 1;

my @dist;

while ($curtime < $time) {
  my $line=<IN>;
  @fields = split(' ', $line);
  $curtime = shift(@fields);
}

my $i = 0;
while (scalar @fields > 0) {
  my $S_i = shift(@fields);
  my $I_i = shift(@fields);
  my $R_i = shift(@fields);
  my $N_i = $S_i + $I_i + $R_i;
  print "$i $N_i\n";
  $i++;
}
