#!/usr/bin/perl

for (my $ninit = 5000; $ninit < 10000; $ninit=$ninit+1000) {
  for (my $beta=0.01;$beta<=0.6;$beta=$beta+0.01) {
    for (my $rho=0.01;$rho<=0.3;$rho=$rho+0.01) {
      for (my $epsilon=1;$epsilon<=5;$epsilon=$epsilon+1) {
        system "./solve_ode.x --file-id angola4 --ic-file angola_$ninit.ic --nsave 1 --dt=0.1 --N 1250000 -p angola3.prm -m angola3.mpm  --beta $beta --rho $rho --sigma-file angola.csv --epsilon $epsilon > /dev/null 2>&1";
        system "./reported.sh angola4.dat > reported2.dat";
        print "$ninit $beta $rho $epsilon ".`./least_squares.sh reported2.dat recorded.dat`;
      }
    }
  }
}
