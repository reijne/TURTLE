#!/bin/sh

aclocal --output=aclocal.tmp
cat aclocal.tmp | sed -e "s/| sed -e 's\/\^lib\/cyg\/'//" > aclocal.tmp.1
mv -f aclocal.tmp.1 aclocal.tmp 
cat m4/* aclocal.tmp > aclocal.m4
rm -f aclocal.tmp
automake
autoconf
