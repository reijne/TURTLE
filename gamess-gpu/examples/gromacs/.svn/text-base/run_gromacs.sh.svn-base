#!/bin/bash
if ( test ${#GROMACS_DIR} -eq 0 ) then
  echo "Please set environment variable GROMACS_DIR to the directory of your"
  echo "GROMACS installation. This test case expects the following executables:"
  echo "  - GROMACS_DIR/src/kernel/grompp"
  echo "  - GROMACS_DIR/src/kernel/mdrun"
  exit 100
fi
mkdir LTvacuo
cd LTvacuo
../create_tops.scr 0.12 0.5 1
../run.scr 0.12 0.5 1
