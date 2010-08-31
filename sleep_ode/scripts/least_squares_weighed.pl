#!/usr//bin/perl

$firstfile = shift(@ARGV);
$secondfile = shift(@ARGV);

(($firstfile ne "") && ($secondfile ne "")) || die "Need to specify two filenames";

my $lq = 0;

my @indices;
my @values;

open (FIRSTFILE, "<$firstfile");

while ($line = <FIRSTFILE>) {
  my @fields = split(/ /,$line);
  push(@indices, $fields[0]);
  push(@values, $fields[4]);
}

close(FIRSTFILE);

open (SECONDFILE, "<$secondfile"); 
my $current_index = 0;
while ($line = <SECONDFILE>) {
  my @fields = split(/\t/,$line);
  if ($fields[0] == $indices[$current_index]) {
#    $lq += sqrt($fields[1]) * ($fields[1] - $values[$current_index])**2;
    $lq += ($fields[1] - $values[$current_index])**2 / $fields[1];
    ++$current_index;
  }
}
close (SECONDFILE);

print "$lq\n";
