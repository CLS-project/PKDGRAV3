#!/bin/sh
#
# This must be done once:
# touch NEWS README AUTHORS ChangeLog
# This should be done when Makefile.am or configure.ac change.  Best to do
# a "make distclean" if already configured before running this script.
test -d openpa || tar -xf openpa.tar.bz2
( cd openpa ; autoheader ; aclocal -Iconfdb ; automake -acf ; autoconf -i -f )
autoheader
aclocal -Im4 -I$(gsl-config --prefix)/share/aclocal
automake -acf
autoconf -i -f
# ./configure
# make
