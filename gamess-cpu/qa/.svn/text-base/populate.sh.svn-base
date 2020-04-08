#!/bin/sh

QADIR=../../GAMESS-UK_misc/QA

if test ! -d  $QADIR
then
  echo You must ensure that the directory $QADIR exists first
  exit 1
fi

(cd $QADIR; ./make_new_qa tar)

echo Extracting QA scripts and files 

tar xf $QADIR/qa.tar

if test ! -f  ./testdirs.txt
then
  echo Now run \"make qa\" in the GAMESS-UK/m4 directory
fi

echo 0
