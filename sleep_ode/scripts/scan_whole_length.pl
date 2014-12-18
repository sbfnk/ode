#!/usr/bin/perl

for (my $ninit = 100; $ninit < 1000; $ninit=$ninit+100) {
  for (my $beta = 1; $beta < 2; $beta = $beta + 0.02) {
    for (my $rho = 0; $rho < 1; $rho = $rho + 0.02) {
      for (my $epsilon = 1; $epsilon <= 5; $epsilon = $epsilon + 1) {
        system "./solve_ode.x --file-id angola4 --ic-file angola_$ninit.ic --nsave 1 --dt=0.1 --N 1250000 -p angola4.prm -m angola4.mpm  --beta $beta --rho $rho --sigma-file angola.csv --epsilon $epsilon > /dev/null 2>&1";
        system "./reported.sh angola4.dat > reported.dat";
        print "$ninit $beta $rho $epsilon ".`./least_squares.sh reported.dat recorded.dat`;
      }
    }
  }
}
#for (my $ninit = 1000; $ninit < 10000; $ninit=$ninit+1000) {
#  for (my $beta = 0; $beta < 1; $beta = $beta + 0.02) {
#    for (my $rho = 0; $rho < 1; $rho = $rho + 0.02) {
#      for (my $epsilon = 1; $epsilon <= 5; $epsilon = $epsilon + 1) {
#        system "./solve_ode.x --file-id angola4 --ic-file angola_$ninit.ic --nsave 1 --dt=0.1 --N 1250000 -p angola4.prm -m angola4.mpm  --beta $beta --rho $rho --sigma-file angola.csv --epsilon $epsilon > /dev/null 2>&1";
#        system "./reported.sh angola4.dat > reported.dat";
#        print "$ninit $beta $rho $epsilon ".`./least_squares.sh reported.dat recorded.dat`;
#      }
#    }
#  }
#}
#for (my $ninit = 10000; $ninit < 100000; $ninit=$ninit+10000) {
#  for (my $beta = 0; $beta < 1; $beta = $beta + 0.02) {
#    for (my $rho = 0; $rho < 1; $rho = $rho + 0.02) {
#      for (my $epsilon = 1; $epsilon <= 5; $epsilon = $epsilon + 1) {
#        system "./solve_ode.x --file-id angola4 --ic-file angola_$ninit.ic --nsave 1 --dt=0.1 --N 1250000 -p angola4.prm -m angola4.mpm  --beta $beta --rho $rho --sigma-file angola.csv --epsilon $epsilon > /dev/null 2>&1";
#        system "./reported.sh angola4.dat > reported.dat";
#        print "$ninit $beta $rho $epsilon ".`./least_squares.sh reported.dat recorded.dat`;
#      }
#    }
#  }
#}
