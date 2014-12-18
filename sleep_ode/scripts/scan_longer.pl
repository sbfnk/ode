#!/usr/bin/perl

for (my $ninit = 30000; $ninit <= 100000; $ninit=$ninit+10000) {
  for (my $beta=0.01;$beta<=0.3;$beta=$beta+0.01) {
    for (my $rho=0.01;$rho<=0.2;$rho=$rho+0.01) {
      for (my $epsilon=1;$epsilon<=5;$epsilon=$epsilon+1) {
        system "./solve_ode.x --file-id angola3 --ic-file angola_$ninit.ic --nsave 1 --dt=0.1 --N 1250000 -p angola3.prm -m angola3.mpm  --beta $beta --rho $rho --sigma-file angola.csv --epsilon $epsilon > /dev/null 2>&1";
        system "./reported.sh angola3.dat > reported.dat";
        print "$ninit $beta $rho $epsilon ".`./least_squares.sh reported.dat recorded2.dat`;
      }
    }
  }
}
