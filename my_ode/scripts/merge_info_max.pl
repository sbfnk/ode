#!/usr/bin/perl

use Getopt::Long;

my $usage = "Usage: merge_finals.pl filename\n";

my $filename = shift(@ARGV);
($filename eq "") && ($usage);


open(IN, "<", $filename) or die "Can't open $filename for reading\n";

my $param = "[0-9].[0-9][0-9]";

my $pattern="./ab_($param)\_($param)/dgl_($param)\_($param)\_($param)/or_($param)\_($param).final\n";
my @data;

while ($dataline = <IN>) {
  $fileline = <IN>;
  if ($fileline =~ /$pattern/) {
    chomp $dataline;
    my @data = split(/ /,$dataline);
    if ($5 > 0) {
      printf "%8.2f", ($data[4]*$6/$5*(1-$6/($6+$4+$5)));
      print " a=$1 b=$2 d=$3 g=$4 l=$5 o=$6 r=$7 $dataline\n";
    }
  }
}
