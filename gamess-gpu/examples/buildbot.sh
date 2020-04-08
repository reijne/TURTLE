#!/bin/bash

# Run the GAMESS-UK tests under buildbot
#
# If the -a argument is supplied, then run all the tests with run_em_all, else just
# run the chap2 examples
#
# if the -p argument is supplied we run the parallel tests and create a file called gamess 
# in the bin directory that invokes an mpirun that is tailored to the machine we are on.
#
# Depending on the machine, this may or may not submit the job to the queue system
#

# For debugging: -v print shell input lines, -x print commands and arguments
#set -x
#set -v

runall=0
parallel=0
# Get command-line options
while getopts "ap" flag $@
do
  case $flag in
      a) runall=1;;
      p) parallel=1;;
      # ? means an invalid option
      ?) echo "$0: Problem with option: -$OPTARG"; exit 1;;
      # : means an option without an argument - aparently not...
      #:) echo "$0: Option $OPTARG requires an argument!" | tee $logfile; exit 1;;
  esac
done

# See what we're running
if [ $runall -eq 1 ]
    then
    if [ $parallel -eq 1 ]
	then
	cmd="./run_em_all.sh -s -p"
    else
	cmd="./run_em_all.sh -s"
    fi
else
    if [ $parallel -eq 1 ]
	then
	cmd="cd ./chap2 && ./run_chap2.sh -s -p"
    else
	cmd="cd ./chap2 && ./run_chap2.sh -s"
    fi
fi

# Source the buildbot utility functions
source ./buildbot_funcs.sh

# Set the location of the binary directory
bindir=`(cd ../bin;pwd)`

#
# Work out where we are and what to run
#
##############################################################################################
#
# cselnx1
#
##############################################################################################
#
if [ `hostname` = "cselnx1" ]
    then
    if [ $parallel -eq 1 ]
	then
	echo '#!/bin/bash' > $bindir/gamess
        #echo mpiexec.gforker -np 4 ${bindir}/gamess-uk '$*' >> $bindir/gamess
        #mpiexec.gforker causes the code to die (why?)
        #using mpirun gets to code to work but there is still a lot of noise
        #complaining about:
        #WARNING: Failed to open "OpenIB-mlx4_0-2" [DAT_PROVIDER_NOT_FOUND:DAT_NAME_NOT_REGISTERED].
        #[0,1,3]: uDAPL on host cselnx1 was unable to find any NICs.
        #the output also contains a lot of:
        #DAT: library load failure: /usr/lib64/libdaplcma.so.1: undefined symbol: dat_registry_add_provider
        #so more stuff to sort out
	echo mpirun -np 4 ${bindir}/gamess-uk '$*' >> $bindir/gamess
	chmod +x $bindir/gamess
     fi
    #
    # Run the tests
    eval $cmd
#
##############################################################################################
#
# csesgi1
#
##############################################################################################
#
elif [ `hostname` = "csesgi1" -a $parallel -eq 1 ]
    then
    echo '#!/bin/bash' > $bindir/gamess
    echo mpirun -np 4 ${bindir}/gamess-uk '$*' >> $bindir/gamess
    chmod +x $bindir/gamess
    eval $cmd
#
#
##############################################################################################
#
# hapu
#
##############################################################################################
#
elif [ `hostname` = "hapu33" -a $parallel -eq 1 ]
    then
    #
    # Create gamess script to run the binary
    #
    echo '#!/bin/bash' > $bindir/gamess
    echo "mpirun -srun $bindir/gamess-uk $\*" >> $bindir/gamess
    chmod +x $bindir/gamess
    #
    # Create submission script for the queueing system
    # This is where to set the number of nodes to run on
    #
    cat > submit.sh <<EOF
#!/bin/sh
# BSUB -n 16
# BSUB -J gamess-uk
# BSUB -o buildbot.log
# BSUB -W 04:00
# BSUB -ext "SLURM[nodes=4]"

./run_em_all.sh -s -p
EOF
    #
    # Set up and load the relevant modules
    #
    . /etc/profile.d/modules.bash
    module load mpi/hp/default
    #
    # Now run using the batch function
    #
    run_lsf submit.sh || exit 1
    #
    # If there is a file called error cat the log to stoud
    [ -f ./ERROR ] && cat ./buildbot.log && exit 1 || exit 0
    #
#
##############################################################################################
#
# plogin
#
##############################################################################################
#
elif [ `hostname` = "plogin" ]
    then
    #
    # Create submission script for the queueing system
    # This is where to set the number of nodes to run on
    #
    cat > submit.sh <<EOF
#!/bin/sh
#@ error   = buildbot.log
#@ output  = buildbot.log
#@ notification  = never
#@ wall_clock_limit=12:00:00
#@ class=par32_12
#@ job_type = serial
#@ node_usage = shared
#@ queue

$cmd
EOF
    #
    # Now run using the batch function
    #
    run_ll submit.sh
    #
    # If there is a file called error cat the log to stoud
    [ -f ./ERROR -o $? -ne 0 ] && cat ./buildbot.log && exit 1 || exit 0
#
##############################################################################################
#
# Just run locally
#
##############################################################################################
#
else
    eval $cmd
fi




#
# HERE BE DRAGONS!
# Old stuff below here
#

# #
# ##############################################################################################
# #
# # HPCx
# #
# ##############################################################################################
# #
# elif [ `hostname` = "l1f401" ]
#     then
#     #
#     # Create gamess script to run the binary
#     #
#     echo '#!/bin/bash' > $bindir/gamess
#     echo "cat > datain" >> $bindir/gamess
#     echo "/usr/bin/poe $bindir/gamess-uk \$*" >> $bindir/gamess
#     chmod +x $bindir/gamess
#     #
#     # Create submission script for the queueing system
#     # This is where to set the number of nodes to run on
#     #
#     cat > submit.sh <<EOF
# #!/bin/sh
# #@ shell = /bin/sh
# #@ job_type = parallel
# #@ account_no = c01-chem
# #@ job_name = VALIDATE
# #@ error = buildbot.log
# #@ output = buildbot.log
# #@ wall_clock_limit = 00:45:00
# #@ tasks_per_node = 16
# #@ cpus = 16
# #@ network.MPI_LAPI = csss,not_shared,IP
# #@ node_usage=not_shared
# #@ queue

# export LAPI_USE_SHM=yes
# export MP_SHARED_MEMORY=yes
# export MP_CSS_INTERRUPT=yes
# export AIXTHREAD_SCOPE=S
# export MP_POLLING_INTERVAL=25000
# export RT_GRQ=ON
# #
# export MP_PULSE=0
# export MP_INTRDELAY=100
# export MP_SINGLE_THREAD=yes
# export VT_ROOT=/usr/local/packages/vampir

# ./run_em_all.sh -s -p
# EOF
#     #
#     # Now run using the batch function
#     #
#     run_ll submit.sh || exit 1
#     #
#     # If there is a file called error cat the log to stoud
#     [ -f ./ERROR ] && cat ./buildbot.log && exit 1 || exit 0
#     #



# #
# ##############################################################################################
# #
# # cseem64t
# #
# ##############################################################################################
# #
# elif [ `hostname` = "cseem64t" ]
#     then
#     #
#     # Create gamess script to run the binary
#     #
#     echo '#!/bin/bash' > $bindir/gamess
#     echo "cat > $bindir/datain" >> $bindir/gamess
#     echo 'PPN=4' >> $bindir/gamess
#     echo "/usr/bin/mpirun -np \$((NSLOTS*PPN)) \
# -machinefile $HOME/.mpich/ndfile.\$JOB_ID  \
# -stdin $bindir/datain \
# $bindir/gamess-uk \$*" >> $bindir/gamess
#     chmod +x $bindir/gamess
#     #
#     # Create submission script for the queueing system
#     # This is where to set the number of nodes to run on
#     #
#     cat > submit.sh <<EOF
# #!/bin/bash
# #$ -pe mpich 2
# #$ -cwd
# #$ -V
# #$ -j yes
# #$ -o run_em_all.log
# #
# # Create the machinefile
# sed s/$/:4/ \$HOME/.mpich/mpich_hosts.\$JOB_ID > \$HOME/.mpich/ndfile.\$JOB_ID
# #
# # Wait 5 minutes to ensure that the operating system has cleaned up the nodes
# # before trying to run anything. Hopefully this will get rid of spurious 
# # crashes due to: "No free InfiniPath contexts available on /dev/ipath".
# #
# sleep 300
# #
# # Run the job
# #
# ./run_em_all.sh -s -p
# EOF
#     #
#     # Now run using the batch function
#     #
#     run_sge submit.sh
#     #
#     # If there is a file called error cat the log to stoud
#     [ -f ./ERROR ] && cat ./run_em_all.log && exit 1 || exit 0
