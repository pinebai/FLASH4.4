!!****if* source/Grid/GridParticles/gr_ptSetIndices
!!
!! NAME
!!
!!  gr_ptSetIndices
!!
!! SYNOPSIS
!!
!!  gr_ptSetIndices(integer(IN) :: index_list(:),
!!                  integer(IN) :: count)
!!
!! DESCRIPTION
!!
!!  GridParticles subunit is used for both, tracking the Lagrangian particles, 
!!  and for tracing rays. The attributes for both data structures, though similar,
!!  do differ some. So the routine within GridParticles cannot make assumptions about
!!  knowing the indices such as positions, tag etc. This routine should be called by 
!!  Particles or RayTrace units to setup the indices correctly so that GridParticles 
!!  doesn't have to know them beforehand.
!!  
!!
!! ARGUMENTS
!!
!! index_list - the list of indices for attributes within the data structures that
!!              need to be known by the GridParticles subunit
!! count      - the count of entries in the list carried by index_list
!!
!!***

#include "GridParticles.h"

subroutine gr_ptSetIndices(index_list,count)

  use gr_ptData, ONLY : gr_ptBlk, gr_ptProc, gr_ptTag, &
       gr_ptPosx, gr_ptPosy, gr_ptPosz,&
       gr_ptPos2x, gr_ptPos2y, gr_ptPos2z,&
       gr_ptVelx, gr_ptVely, gr_ptVelz,&
       gr_ptVel2x, gr_ptVel2y, gr_ptVel2z,&
       gr_ptPosTmp, gr_ptVel, gr_ptVelTmp,gr_ptVirtual

  implicit none
  integer, intent(IN) :: count
  integer,dimension(count), intent(IN) :: index_list


  gr_ptPosx = index_list(GRPT_POSX_IND)
  gr_ptPosy = index_list(GRPT_POSY_IND)
  gr_ptPosz = index_list(GRPT_POSZ_IND)
  gr_ptPos2x = index_list(GRPT_POSXTMP_IND)
  gr_ptPos2y = index_list(GRPT_POSYTMP_IND)
  gr_ptPos2z = index_list(GRPT_POSZTMP_IND)
  gr_ptVelx = index_list(GRPT_VELX_IND)
  gr_ptVely = index_list(GRPT_VELY_IND)
  gr_ptVelz = index_list(GRPT_VELZ_IND)
  gr_ptVel2x = index_list(GRPT_VELXTMP_IND)
  gr_ptVel2y = index_list(GRPT_VELYTMP_IND)
  gr_ptVel2z = index_list(GRPT_VELZTMP_IND)
  gr_ptBlk  = index_list(GRPT_BLK_IND)
  gr_ptProc = index_list(GRPT_PROC_IND)
  gr_ptTag  = index_list(GRPT_TAG_IND)
  gr_ptPosTmp = (index_list(GRPT_POSTMP)==GRPT_EXIST)
  gr_ptVel = (index_list(GRPT_VEL)==GRPT_EXIST)
  gr_ptVelTmp = (index_list(GRPT_VELTMP)==GRPT_EXIST)
  gr_ptVirtual = index_list(GRPT_VIRTUAL)

end subroutine gr_ptSetIndices
