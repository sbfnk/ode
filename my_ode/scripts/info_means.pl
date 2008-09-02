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
my $sum = .0;
my $sqsum = .0;
my $count = 0;

while (scalar @fields > 0) {
  my $S_i = shift(@fields);
  my $I_i = shift(@fields);
  my $R_i = shift(@fields);
  my $N_i = $S_i + $I_i + $R_i;
  if (scalar @fields > 0) {
    $sum += $N_i * $i;
    $sqsum += $N_i * $i*$i;
    $count += $N_i;
  }
  $i++;
}

my $mean;
my $dev;

if ($count > 0) {
  $mean = $sum / $count;
  $dev = sqrt(($sqsum - 2*$sum*$mean + $count*$mean*$mean) / $count);
} else {
  $mean = 0;
  $dev = 0;
}

print "$mean $dev\n";

