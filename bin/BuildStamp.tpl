## Lines starting with ## are comments inside template file
## All other lines including empty lines are non-comments
## 
## This file is a template for generating the make_bstamp script.
## For syntax of this file see "Readme.template".
##
## VALID VARIABLE NAMES FOR THIS TEMPLATE
##
## date  -> date/time when setup script was run
## uname -> Machine information tuple (usually length 4)
##
##
#!/bin/sh -f

#\t\tScript to create a subroutine which writes a 'build stamp'
#\t\tto the log file

rm -f setup_buildstamp.F90
echo '      ' subroutine setup_buildstamp \(s_stamp_str, b_stamp_str, str_len\) >> setup_buildstamp.F90
echo '      ' implicit none >> setup_buildstamp.F90
echo '      ' integer                  :: str_len >> setup_buildstamp.F90
echo '      ' character\(len=str_len\) :: s_stamp_str, b_stamp_str >> setup_buildstamp.F90
echo '      ' s_stamp_str = \\'%(date)s\\'  >> setup_buildstamp.F90
echo '      ' b_stamp_str = \\'`date '+%%a %%b %%e %%H:%%M:%%S %%Y'`\\'  >> setup_buildstamp.F90
echo '      ' return >> setup_buildstamp.F90
echo '      ' end subroutine >> setup_buildstamp.F90
echo '      ' >> setup_buildstamp.F90
echo '      ' subroutine setup_systemInfo \(system_str, str_len\) >> setup_buildstamp.F90
echo '      ' integer                  :: str_len >> setup_buildstamp.F90
echo '      ' character\(len=str_len\) :: system_str >> setup_buildstamp.F90
echo '      ' system_str = "'%(uname!&\n& )s'"  >> setup_buildstamp.F90
echo '      ' return >> setup_buildstamp.F90
echo '      ' end subroutine >> setup_buildstamp.F90
echo '      ' >> setup_buildstamp.F90

