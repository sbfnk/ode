#!/usr/bin/perl

my $codeDir = $ENV{'CODEDIR'}."/sleep_ode";
my $dataDir = $ENV{'DATADIR'}."/Sleeping_Sickness";

for (my $nmul = 1; $nmul < 4; $nmul=$nmul+1) {
  for (my $ninit = 100 * (10**($nmul - 1)); $ninit < 1000 * (10**($nmul - 1)); $ninit=$ninit+(100 * (10**($nmul - 1)))) {
    for (my $beta=(0.5*(2**(3-$nmul)));$beta<(1*(2**(3-$nmul)));$beta=$beta+0.02) {
#      for (my $rho=0.02;$rho<=1;$rho=$rho+0.02) {
      for (my $rho=0.02;$rho<2;$rho=$rho+0.02) {
        for (my $epsilon=0.1;$epsilon<=10;$epsilon=$epsilon*10) {
          system "$codeDir/bin/solve_ode.x --file-id $dataDir/temp --ic-file $codeDir/init/angola_$ninit.ic --nsave 1 --dt=0.1 --N 1250000 -p $codeDir/params/angola.prm -m $codeDir/params/angola2.mpm  --t0=1972 --tmax=1998 --beta $beta --rho $rho --sigma-file $codeDir/data/screening.dat --epsilon $epsilon > /dev/null 2>&1";
          system "$codeDir/scripts/reported.sh $dataDir/temp.dat > $dataDir/temp_reported.dat";
	  my $distance = `$codeDir/scripts/least_squares_weighed.pl $dataDir/temp_reported.dat $codeDir/data/recorded.dat`;
	  chomp $distance;
          my $infected_total_last = 
            `tail -n 1 $dataDir/temp.dat | awk '{print \$3" "\$5}'`;
          chomp $infected_total_last;
          print "$ninit $beta $rho $epsilon $infected_total_last $distance\n";
        }
      }
    }
  }
}
