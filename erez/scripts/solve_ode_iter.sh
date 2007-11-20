#!/bin/bash

# usage message
usage="Usage: odesolve.sh file_id [-m model_params_file] [-o ode_params_file]"
usage="$usage [-ic ic_file] [-f] [-v] [-c] [-plot] [-gv] [-mf/-pa]"
usage="$usage [other_options, e.g. tmax, dt, etc.]"

# saving command line 
comm_line="$0 $@"

# set file_id
file_id=$1
shift

# check that file_id is not NULL
if [ -z $file_id ]; then
    echo "ERROR: file_id not set"
    echo "$usage"
    exit 1
fi

# set path
ode_dir=$CODEDIR/ode

# init message
echo
echo 'Running ode solver script'
echo '========================='
echo

# defaults
echo '... setting defaults'
force=n
verbose=n
convergence=n
plot=n
gv=n
model_type='mf'

# parse argument list
echo '... parsing args'
while [ $# -ge 1 ]; 
  do
  
  case $1 
      in
      -m) shift; model_prm=$1;;
      -o) shift; ode_prm=$1;;
      -ic) shift; ic_file=$1;;
      -f) force="y";;
      -v) verbose="y";;
      -c) convergence="y";;
      -plot) plot="y";;
      -h) echo $usage; exit 0;;
      -gv) gv="y";;
      -pa) model_type='pa';;
      -mf) model_type='mf';;
      *) ode_options="$ode_options $1";;
  esac
  shift
  
done

# reset variables and path
echo '... resetting variables and paths'
solver=solve_"$model_type"_ode.x 
exec=$ode_dir/bin/$solver

file_id="$file_id"_$model_type
output_dir="$DATADIR/$file_id"
log_file=$file_id.log

# check if output_dir exists
echo '... checking output_dir'
if [ -d $output_dir ]; then
    if [ $force == "y" ]; then
	echo "... overwriting contents of $output_dir"
        rm -rf $output_dir/*
    else
	echo "ERROR: $output_dir exists"
	echo "Use -f to override"
	exit 1
    fi
else
    mkdir $output_dir
fi

# writing comm line
echo $comm_line > $output_dir/$file_id.comm_line

# check for files before copy
echo "... copy files to output_dir"

# exec file
cp -f $exec $output_dir
exec=$output_dir/$solver

# ic_file
if [ -e $ic_file ]; then
    cp -f $ic_file $output_dir/$file_id.init
else
    echo "ERROR: ic_file $ic_file does not exist"
    exit 1
fi

# ode_prm
if [ -e $ode_prm ]; then
    cp -f $ode_prm $output_dir/$file_id.ode.prm
else
    echo "ERROR: ode_prm $ode_prm does not exist"
    exit 1
fi

# model_prm
if [ -e $model_prm ]; then
    cp -f $model_prm $output_dir/$file_id.model.prm
else
    echo "ERROR: model_prm $model_prm does not exist"
    exit 1
fi

# extract_pairs.m
#if [ -e extract_pairs.m ]; then
#    cp -f extract_pairs.m $output_dir/$file_id.m
#fi

# set ode command line options
echo "... setting command line options for ode solver"
options="--file-id=$output_dir/$file_id"
options="$options -o $output_dir/$file_id.ode.prm"
options="$options -m $output_dir/$file_id.model.prm"
options="$options $ode_options"

if [ $verbose == "y" ]; then
    options="$options --verbose"
fi

if [ $convergence == "y" ]; then
    options="$options --check-convergence"
fi

# print message
echo
echo "... Running ode solver for $file_id"
echo "==============================================="
####################################################
# loop 

# setting dat file
dat_file=chi.tau.$model_type.dat

for chi in $(seq 0 .1 10)
  do
  chi_values="--chi=$chi"

  # print message
  #echo "chi = $chi"

  for tau in $(seq 0 .1 5)
    do
    
    # print message
    echo "chi = $chi : tau = $tau"
    
    # rescale tau+- t++ by sigma
    tau_sig=$(echo "scale=8; 0.1 * $tau" | bc)
    
    # setting tau
    tau_values="--tau--=$tau --tau-+=$tau --tau+-=$tau_sig --tau++=$tau_sig"

    # setting ode solver command
    ode_command="$exec $options $tau_values $chi_values"
    
    # execute ode solver
    $ode_command >> $output_dir/$log_file
    echo ' ' >> $output_dir/$log_file
    echo '===============================================================' >> $output_dir/$log_file
    echo ' ' >> $output_dir/$log_file
    
    # sed log file to remove the LONG directory prefix
    tmpdir="$(echo $DATADIR | sed -e 's/\//\\\//g')\/$file_id\/"
    sed -e "s/$tmpdir//g" $output_dir/$log_file > tmp.$file_id
    mv -f tmp.$file_id  $output_dir/$log_file
    
    # extract last line
    data=`cat $output_dir/$file_id.final`; data="$chi $tau $data"
    echo $data >> $output_dir/$dat_file
    
  done
done
####################################################

echo "... done"
echo


exit 0

# plot
  if [ $plot == "y" ]; then
      echo "Making plots"
      
      if [ $nosim == "n" ]; then
	  plot_command="./makeplot.bash $file_id sim"
      else
	  plot_command="./makeplot.bash $file_id nosim"
      fi
      
      $plot_command
      
    # gv
      if [ $gv == "y" ]; then
	  gv $output_dir/$file_id.ps &
      fi
      
  fi
  
exit 0

### example
# ./solve_ode.sh test201 -m model.prm -o ode.prm -ic test34.init -v -c --tmax=500 --vertices 10000 -f --omega=0.01 --nu=0.01 &