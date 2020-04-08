#! /bin/sh

# SWIG is used to generate wrappers for Python and Perl

echo "Generating Python wrapper..."

swig -python -module libpyagentx -Iinclude -o python/pyagentx.tmp python/pytypemaps.i

# dealing with the definition of _GNU_SOURCE

cat python/pyagentx.tmp | sed -e "s/include <Python.h>/include <Python.h>\n#ifdef HAVE_CONFIG_H\n#include<config.h>\n#endif/g" > python/pyagentx.c
rm python/pyagentx.tmp

echo "Generating Perl wrapper..."
swig -perl -module libplagentx -Iinclude -o perl/plagentx.tmp perl/pltypemaps.i

#dealing with the definition of _GNU_SOURCE

echo "#define _GNU_SOURCE" > perl/plagentx.c
cat perl/plagentx.tmp >> perl/plagentx.c
rm perl/plagentx.tmp

echo "Generating Java wrapper..."

swig -java -module libjagentx -Iinclude -o java/jagentx.c java/jtypemaps.i