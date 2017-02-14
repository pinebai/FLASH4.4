#!/bin/sh -f

# script to create a brief summary of svn situation
# to append to FLASH release

do_cdup=""

rawversions=$(echo \
 `svn status -v ../source ../bin|grep -v '^[?]  '|grep -v '^Status '|cut -b 9-|sed -e 's/^ *//'| \
  sed -e '/-/d' -e 's/ .*//'|sort -n -r|sed -e 's/^/r/'|uniq -c|sort -n -r|sed -e '1s/^ *[0-9]* */svn:/'| \
  tr -d ' '`|tr ' ' ',')

if [ -z "$rawversions" ]; then
    DIRUP=$(dirname $(readlink -e $(dirname ../source/..)))
    if [ -n "$DIRUP" ]; then
      do_cdup=y
      rawversions=$(echo \
	`svn status -v "$DIRUP"/source "$DIRUP"/bin|grep -v '^[?]  '|grep -v '^Status '|cut -b 9-|sed -e 's/^ *//'| \
	sed -e '/-/d' -e 's/ .*//'|sort -n -r|sed -e 's/^/r/'|uniq -c|sort -n -r|sed -e '1s/^ *[0-9]* */svn:/'| \
	tr -d ' '`|tr ' ' ',')
    fi
fi

versions=$rawversions

if [ -z "$do_cdup" ]; then
    info=$(echo `svn st ../source ../bin|grep -v '^[?]  '|cut -b 1|sort|uniq -c|tr -d ' '`|tr ' ' ',')
else
    info=$(echo `svn st "${DIRUP}/"source "${DIRUP}/"bin|grep -v '^[?]  '|cut -b 1|sort|uniq -c|tr -d ' '`|tr ' ' ',')
fi
if [ -n "$info" ]; then
    versions="$rawversions,changed:$info"
fi

if [ -n "$versions" ]; then
    echo "($versions)"
fi
