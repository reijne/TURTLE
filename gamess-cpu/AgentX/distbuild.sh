#! /bin/sh

if ! test $# -eq 1 ; then

echo "useage: distbuild 1|0"
exit 1

fi

# echo "Preparing distribution"

./autogen.sh
cd libxml2
./autogen.sh
cd ..

./configure --prefix=`pwd`/local --with-python --with-perl --with-fortran --with-axtransform --with-libxml2
make uninstall
make maintainer-clean

# set variables

PACKAGE_RELEASE=`cat RELEASE`

if test $1 -eq 1 ; then
PACKAGE_NAME="AgentX-libxml2"
PACKAGE_OPTS="--with-libxml2"
PACKAGE_DEPS="glibc >= 2.3"
else
PACKAGE_NAME="AgentX"
PACKAGE_OPTS="--without-libxml2"
PACKAGE_DEPS="glibc >= 2.3, libxml2 >= 2.6"
fi

echo "Preparing license"

cat LICENSE | sed -e "s/RELEASE [0-9].[0-9].[0-9a-z]/RELEASE ${PACKAGE_RELEASE}/g" > LICENSE.dist
mv -f LICENSE.dist LICENSE

echo "Preparing RPM spec file"

cat agentx.spec.in | sed -e "s/@PACKAGE_NAME@/${PACKAGE_NAME}/g" -e "s/@PACKAGE_RELEASE@/${PACKAGE_RELEASE}/g" -e "s/@PACKAGE_OPTS@/${PACKAGE_OPTS}/g" -e "s/@PACKAGE_DEPS@/${PACKAGE_DEPS}/g" > agentx.spec

echo "Preparing tar archive"

cd ..
tar czvf AgentX-${PACKAGE_RELEASE}.tar.gz --exclude .svn AgentX-${PACKAGE_RELEASE}

echo "Building RPM"

rpmbuild -tb AgentX-${PACKAGE_RELEASE}.tar.gz

exit 0
