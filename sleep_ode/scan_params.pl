#!/usr/bin/perl

for (my $beta=0;$beta<1;$beta=$beta+0.01) {
  for (my $rho=0;$rho<1;$rho=$rho+0.01) {
    system "./solve_ode.x --file-id angola3 --ic-file angola3.ic --nsave 1 --dt=0.1 --N 2500000 -p angola3.prm -m angola3.mpm  --beta $beta --rho $rho > /dev/null 2>&1";
    system "./reported.sh angola3.dat > reported.dat";
    print "$beta $rho ".`./least_squares.sh reported.dat recorded.dat`;
  }
}
