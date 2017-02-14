#ifndef _GR_PTMAPTOMESH_H
#define _GR_PTMAPTOMESH_H

#include "constants.h"

#if 0
*************************************************************************
The Paramesh surr_blks data structure is indexed using the following labels.
*************************************************************************
#endif

#define BLKID BLKNO
#define BLKPROC PROCNO
#define BLKTYPE TYPENO



#if 0
*************************************************************************
The negh() data structure is indexed using the BLKID, BLKPROC, and 
  REFLEVELDIF labels.
  *************************************************************************
#endif

#define REFLEVELDIF TYPENO
#define SIZENEGH REFLEVELDIF



#if 0
  *************************************************************************
  The metadata header is indexed using BLKID, CORNERID and COORDSID.
  e.g.
  3D simulation: BLKID=1, CORNERID=4, COORDSID=7, SIZE_HEADER=12
  *************************************************************************
#endif

#define CORNERID (REFLEVELDIF + 1)
#define COORDSID (CORNERID + MDIM)
#define SIZE_HEADER (COORDSID + (MDIM*2) - 1)


#if 0
  *************************************************************************
  ABSMAXNEGH is the same as 2**(MDIM-1).  It is required to prevent 
  compile time errors in the routine gr_ptFindNegh.  Here, there are  
  compile time problems when MAXNEGH is used instead.
  *************************************************************************
#endif

#ifndef ABSMAXNEGH
#define ABSMAXNEGH 4
#endif

#endif
