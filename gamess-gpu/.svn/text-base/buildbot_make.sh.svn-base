#!/bin/bash

# Build GAMESS-UK under buildbot
#
# This just ensure that the environment is correct to build the code on the specified machine
#
#
##############################################################################################
#
# hapu
#
##############################################################################################
#
if [ `hostname` = "hapu33" ]
    then
    . /etc/profile.d/modules.bash
    module load ifort/11.0/default icc/11.0/default mpi/hp/default
fi

# Assume the environment now o.k. so run make
make
