#!/usr/bin/sh
cd /tmp/edl/
doxygen
cp ./logo.png /srv/www/htdocs/edl
cd /tmp
rm -rf edl
