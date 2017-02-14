#!/usr/bin/env bash

# Generic file handling routine
# Takes filename as only argument
function handle()
{
  FILE=$1
  # create the .orig if need be (so .orig is the link to repos)
  if [ ! -s $FILE.orig ]; then
     mv $FILE $FILE.orig;
  fi
  # if orig has changed or $FILE is not there, regenerate $FILE 
  if [ \( ! -e $FILE \) -o \( $FILE.orig -nt $FILE \) ]; then
     ./setup_reorder.py --input=$FILE.orig --output=$FILE --auto;
  fi
}

# special handling for Paramesh
function handle_PM()
{
  FILE=$1
  # create the .orig if need be (so .orig is the link to repos)
  if [ ! -s $FILE.orig ]; then
     mv $FILE $FILE.orig;
  fi
  # if orig has changed or $FILE is not there, regenerate $FILE 
  if [ \( ! -e $FILE \) -o \( $FILE.orig -nt $FILE \) ]; then
     ./setup_reorder.py --input=$FILE.orig --output=$FILE --auto --five=unk1 --five=tgvar[xyz] --five=t_unk --five=gt_unk --four=recv --four=unk1_fl --five=flux_[xyz] --five=tflux_[xyz] --five=ttflux_[xyz] --four=bndtemp[xyz]1 --five=gt_facevar[xyz] --four=ttunk --four=ttfacevarx --four=ttfacevary --four=ttfacevarz --five=facevar[xyz]1 --four=facevar[xyz]1_fl --four=recvf --four=sendf --four=tempf --four=recv[xyz] --four=recv[xyz]f --five=unk_e_[xyz] --five=unk_e_[xyz]1 --five=gt_unk_e_[xyz] --five=unk_n --five=unk_n1 --five=gt_unk_n --five=bedge_facex_[yz] --five=bedge_facey_[xz] --five=bedge_facez_[xy] --five=tbedge_facex_[yz] --five=tbedge_facey_[xz] --five=tbedge_facez_[xy] --four=unk_e_[xyz]1_fl --four=recvar[xyz]1e --four=recvar[xyz]2e --four=recve --four=recvn --four=recvf[xyz] --four=ttunk_e_[xyz] --four=ttunk_n --five=t_unk_e_[xyz] --five=t_unk_n --four=unk_n1_fl --four=temp --four=recvn0 --four=tempn --four=sendn --four=datain --four=dataout --four=unkt --four=facevar[xyz]t --four=unk_e_[xyz]t --four=unk_nt;
  fi
}



function unhandle()
{
  FILE=$1
  # remove $FILE and move .orig to $FILE

  # no original file, nothing to do
  if [ ! -s $FILE.orig ]; then
     return;
  fi

  # if new file exists remove it
  if [ -s $FILE ]; then
     rm -f $FILE 2>/dev/null;
  fi

  # now move the file over
  mv $FILE.orig $FILE
}

# list of files to be processed
# special handling required for Paramesh
FILELIST='%(reordlist)s'
FILELIST_PM='%(reordlist_pm)s'

OPTION=$1

if [ -z $OPTION ]; then
   for m in $FILELIST; do
       handle $m;
   done;
   for m in $FILELIST_PM; do
       handle_PM $m
   done;
elif [ "$OPTION" == "--clean" ]; then
   for m in $FILELIST; do
       unhandle $m;
   done;
   for m in $FILELIST_PM; do
       unhandle $m;
   done;
else
   echo "Unknown option: [$OPTION]"
   exit 1;
fi

