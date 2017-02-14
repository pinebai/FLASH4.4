#!/usr/bin/env bash

START_URL="http://flash.uchicago.edu/website/codesupport/secure/flash3_ug/"
WGET_FLAGS="-r -nv -nH --cut-dirs=1 -p -I website/codesupport" 

# run wget recursively, constrain to codesupport
# filter out all non 404 messages and send the rest
# to std out

wget ${START_URL} ${WGET_FLAGS} 2>&1 | grep -B 1 "404"


