#!/bin/sh
#
# This must be done once:
# touch NEWS README AUTHORS ChangeLog
# This should be done when Makefile.am or configure.ac change.  Best to do
# a "make distclean" if already configured before running this script.
autoheader
aclocal
autoconf -i
automake -a
# ./configure
# make
