#
# Use
#    source THIS_FILE
# or
#    . THIS_FILE
# from bash; other Bourne shells or bash invoked as /bin/sh are probably
# lacking some required features.

pathmungeany () {
# Usage: pathmungeany COMPONENT [WHERE] [VARNAME]
#  COMPONENT :- the component(s) to (maybe) add to the path
#  WHERE :- "after" means: append unless present anywhere
#           "first" means: prepend unless already present at beginning
#           anything else means: prepend unless already present anywhere
#  WHERE :- the name of the environment variable to test and (maybe) modify
#           if empty, defaults to PATH
# Examples:
#    pathmungeany ${MPI}/bin
#    pathmungeany ${MPI}/lib first LD_LIBRARY_PATH
#    pathmungeany ${MPI}/share/man before MANPATH
# Notes:
#  This functions expects that COMPONENT and VARNAME (if given and nonempty)
#  have sane values.  Shit may well happen if special characters like
#  '$', ';', '(', ')', '|', etc. etc. are present.
    if [ -z "$3" ]; then
	varname="PATH"
    else
	varname="$3"
    fi
    if [ -f /bin/egrep -a -x /bin/egrep ]; then
	EGREP=/bin/egrep
    else
	EGREP=/usr/bin/egrep
    fi
    # ${!var}: indirect variable expansion, requires bash-2.0 or later
    if [ -z "${!varname+set}" ] ; then
        eval "export ${varname}"'=$1'
    elif [ "$2" = "first" ] && 
	! echo "${!varname}" |($EGREP -q  "^$1($|:)") ; then
        eval "export ${varname}"'=$1:${'"${varname}"'}'
    elif ! echo "${!varname}" |($EGREP -q "(^|:)$1($|:)") ; then
        if [ "$2" = "after" ] ; then
            eval "export ${varname}"'=${'"${varname}"'}:$1'
        else
            eval "export ${varname}"'=$1:${'"${varname}"'}'
        fi
    fi
}
