#!/usr/bin/perl

my @params=(0.05,0.15,0.25,0.3,0.35,0.4,0.45,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95);
my @omparams=(0.01,0.1,0.2);
my @rhoparams=(0.5,0.8,0.9,0.95);
for ($a=0;$a<(scalar @params);$a=$a+1) {
  $alpha=$params[$a];
  $beta=1;
  $dir1=sprintf("ab_%.2f_%.2f",$alpha,$beta);
  $gamma=0.5;
  $lambda=0.5;
  $dir2=sprintf("gl_%.2f_%.2f",$gamma,$lambda);
  for ($o=0;$o<(scalar @omparams);$o=$o+1) {
    $omega=$omparams[$o];
    for ($r=0;$r<(scalar @rhoparams);$r=$r+1) {
      $rho=$rhoparams[$r];
      $filename=sprintf("~/Data/ODE/SIR/$dir1/$dir2/or_%.2f_%.2f",$omega,$rho);
      $param_part=sprintf("%.2f %.2f %.2f %.2f %.2f %.2f",$alpha,$beta,$gamma,$lambda,$omega,$rho);
      system("echo -ne $param_part \\ ");
      system("./analyse_ode_final.pl $filename\.final | awk '{print \$3}'");
    }
  }
}
