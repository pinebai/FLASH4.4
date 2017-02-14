#!/bin/bash
#
# A simple script which pretty prints the output from the unix
# size utility and then calculates the maximum stack+heap memory
# available to the executable in each BG/P mode.
# It is used as follows:
#   binary_size.sh ./flash4
#
# The script also displays the N biggest bss memory allocations, 
# where N is 10 by default.  This can be changed by
# passing an optional second argument as follows:
#   binary_size.sh ./flash4 20
# 

DEFAULT_NUM_BSS_ALLOCS=10

if [ $# -lt 1 ]; then
  echo "Usage: $0 binary [num_bss_allocs]"
  exit -1
fi

# Store script arguments in bash variables.
binary=$1
if [ $# -eq 1 ]; then
    num_bss_allocs=${DEFAULT_NUM_BSS_ALLOCS}
else
    num_bss_allocs=$2
fi

# Print summary information about the executable.
size ${binary} | awk '{ if (NR==2) \
{ tx=($1/1024^2); da=($2/1024^2); bs=($3/1024^2); sz=($4/1024^2); \
printf "\
Running size on executable:\n\
 text %10.2f MB\n\
 data %10.2f MB\n\
  bss %10.2f MB\n\
 ------------------\n\
 size %10.2f MB\n\n\
This means that the stack + heap on BG/P cannot exceed:\n\
  SMP %10.2f MB\n\
   CO %10.2f MB\n\
   VN %10.2f MB\n\n", \
tx, da, bs, sz, 2048-sz, 1024-sz, 512-sz; } \
}'

# Print the biggest bss memory allocations in the executable.
echo "The ${num_bss_allocs} biggest bss memory allocations in MB"
objdump -t --section=.bss ${binary} | \
    awk '{ if ($4 == ".bss") printf "%-10.6f %s\n", strtonum("0x" $5)/1024^2, $6}' | \
    sort -g -k 1 -r | head -${num_bss_allocs}
