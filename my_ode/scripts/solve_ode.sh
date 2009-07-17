#!/bin/bash

function solve()
{
    # execute ode solver
    $ode_command >> $output_dir/$log_file
    echo ' ' >> $output_dir/$log_file
    echo '===============================================================' >> $output_dir/$log_file
    echo ' ' >> $output_dir/$log_file
    
    # sed log file to remove the LONG directory prefix
    tmpdir="$(echo $DATADIR | sed -e 's/\//\\\//g')\/$file_id\/"
    sed -e "s/$tmpdir//g" $output_dir/$log_file > tmp.$file_id
    mv -f tmp.$file_id  $output_dir/$log_file
    
    return
}

# usage message
usage="Usage: ode_solve.sh file_id [-m model-params-file] [-o ode-params-file]"
usage="$usage [-g graph-params-file] [-ic ic_file] [-f] [-v] [-c]"
usage="$usage [other_options, e.g. tmax, dt, ode-type, etc.]"
usage="$usage [-var1 arg start step stop]"
usage="$usage [-var2 arg start step stop]"


# saving command line 
comm_line="$0 $@"

# set file_id
file_id=$1
shift

# check if -h
if [ $file_id == "-h" ]; then
    echo $usage
    exit 1
fi

# check that file_id is not NULL
if [ -z $file_id ]; then
    echo "ERROR: file_id not set"
    echo "$usage"
    exit 1
fi

# init message
echo
echo 'Running ode solver script'
echo '========================='
echo

# defaults
echo '... setting defaults'

# parse argument list
echo '... parsing args'
while [ $# -ge 1 ]; 
  do
  
  case $1 
      in
      -m) shift; model_prm=$1;;
      -o) shift; ode_prm=$1;;
      -g) shift; graph_prm=$1;;
      -ic) shift; ic_file=$1;;
      -f) force=1;;
      -v) verbose=1;;
      -c) convergence=1;;
      -iter) iter=1;;
      -var1) shift; var1=$1; shift; var1_start=$1; shift; var1_step=$1; shift; var1_stop=$1;;
      -var2) shift; var2=$1; shift; var2_start=$1; shift; var2_step=$1; shift; var2_stop=$1;;
      -h) echo $usage; exit 0;;
      *) ode_options="$ode_options $1";;
  esac
  shift
  
done

# setting variables and paths
echo '... resetting variables and paths'
solver=solve_ode.x 
exec=$CODEDIR/my_ode/bin/$solver

output_dir="$DATADIR/$file_id"
log_file=$file_id.log

# check if output_dir exists
echo '... checking output_dir'
if [ -d $output_dir ]; then
    if [ $force ]; then
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

# graph_prm
if [ -e $graph_prm ]; then
    cp -f $graph_prm $output_dir/$file_id.graph.prm
else
    echo "ERROR: graph_prm $graph_prm does not exist"
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
options="$options -g $output_dir/$file_id.graph.prm"
options="$options --ic-file $output_dir/$file_id.init"
options="$options $ode_options"

if [ $verbose ]; then
    options="$options --verbose"
fi

if [ $convergence ]; then
    options="$options --check-convergence"
fi

# setting ode solver command
ode_command="$exec $options
"
# print message
echo
echo "... Running ode solver for $file_id"
echo "==============================================="
echo
####################################################
# loop 

# print message
if [ $var1 ]; then
    echo -n "Iterating $var1 in [$var1_start:$var1_step:$var1_stop]"
    if [ $var2 ]; then
	echo " and $var2 in [$var2_start:$var2_step:$var2_stop]"
    fi

    if [ $var2 ]; then
	echo -n
    else 
	echo
    fi
fi

if [ $var1 ]; then
    
    # setting dat file
    dat_file=$file_id.$var1.dat
    
    for ivar1 in $(seq $var1_start $var1_step $var1_stop)
      do
      
      # print message
      if [ $var2 ]; then
	  echo -n
      else
	  echo "$var1 = $ivar1"
      fi
      
      # setting ode solver command
      ode_command="$exec $options --$var1=$ivar1"
      
      if [ $var2 ]; then

          # setting dat file
	  dat_file=$file_id.$var1.$var2.dat
	  
	  for ivar2 in $(seq $var2_start $var2_step $var2_stop)
	    do
	    
            # print message
	    echo "$var1 = $ivar1 : $var2 = $ivar2"

            # setting ode solver command
	    if [ $var2 == "tau" ]; then
		ode_command="$exec $options --$var1=$ivar1 --$var2--=$ivar2 --$var2-+=$ivar2"
	    else
		ode_command="$exec $options --$var1=$ivar1 --$var2=$ivar2"
	    fi
	    
            # execute ode solver
	    solve # a function
	    
            # extract last line
	    data=`cat $output_dir/$file_id.final`; data="$ivar1 $ivar2 $data"
	    echo $data >> $output_dir/$dat_file
	    
	  done
	  
      else
	  
          # execute ode solver
	  solve # a function
	  
          # extract last line
	  data=`cat $output_dir/$file_id.final`; data="$ivar1 $data"
	  echo $data >> $output_dir/$dat_file
	  
      fi
      
    done
else    
    solve
fi
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
erez@blin01 ~ $
erez@blin01 ~ $ cd $CODEDIR/solver/scripts
erez@blin01 ~ $ ./solve_ode_iter.sh tmp -mf -m ../params/model.prm -o ../params/ode.prm -g ../../data/graphs/rrg/N.1e4.Qd.6.Qi.5.olp.2/N.1e4.Qd.6.Qi.5.olp.2.graph_stats -ic ../../data/graphs/rrg/N.1e4.Qd.6.Qi.5.olp.2/N.1e4.Qd.6.Qi.5.olp.2.init -v -c -f --tmax=40 --omega=0.03 --mu=0.04 --nsave=10 

detailed:
./solve_ode_iter.sh tmp 
-m ../params/model.prm 
-o ../params/ode.prm 
-g ../../data/graphs/rrg/N.1e4.Qd.6.Qi.5.olp.2/N.1e4.Qd.6.Qi.5.olp.2.graph_stats 
-ic ../../data/graphs/rrg/N.1e4.Qd.6.Qi.5.olp.2/N.1e4.Qd.6.Qi.5.olp.2.init 
-v -c -f 
--tmax=40 --omega=0.03 --mu=0.04 --nsave=10 
-var1 chi 1.2 0.1 1.4 (optional)
-var2 mu 2.1 0.01 2.13 (optional)
