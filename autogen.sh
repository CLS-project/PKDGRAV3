#!/bin/sh
#
autoheader
# touch NEWS README AUTHORS ChangeLog
# touch stamp-h
aclocal
autoconf -i
automake -a
# ./configure
# make
