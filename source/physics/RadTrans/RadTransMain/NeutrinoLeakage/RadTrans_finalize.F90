!!****if* source/physics/RadTrans/RadTransMain/NeutrinoLeakage/RadTrans_finalize
!!
!! NAME
!!
!!  RadTrans_finalize
!!
!! SYNOPSIS
!!
!!  call RadTrans_finalize ()
!!
!! DESCRIPTION
!!
!!  Cleans up the RadTrans unit.
!!
!! ARGUMENTS
!!
!!  NOTES
!!      This unit implements ray-by-ray multispecies neutrino leakage.
!!      Parts of this unit are released under a different license than the
!!      usual FLASH license.  Specifically, some subroutines in rt_calcLeak.F90 and 
!!      rt_calcTau.F90 are released under the Creative Commons 
!!      attribution-noncommercial-share alike license.  Basically, if you use this
!!      unit in your work, the license requires that you cite the two articles 
!!      mentioned below.  More details may be found here:  stellarcollapse.org/codes.html.
!!
!!      * O'Connor, E.P., & Ott, C.D. 2010, CQGra, 27, 114103
!!      * Couch, S.M., & O'Connor, E.P. 2013, arXiv:1310.5728
!!
!!***
subroutine RadTrans_finalize ()
#include "Flash.h"

  use rt_data
  use RadTrans_data, ONLY : rt_useRadTrans

  implicit none
  
  if (.NOT. rt_useRadTrans) return

#ifndef LEAK_STATIC
  if (Allocated(rt_leakArr)) &
       deallocate(rt_leakArr)
  if (Allocated(rt_leakSrc)) &
       deallocate(rt_leakSrc)

  deallocate(rt_leakRadii)
  deallocate(rt_leakTheta)
  deallocate(rt_leakX)
  deallocate(rt_leakY)

  deallocate(rt_dr)
#endif

  return
end subroutine RadTrans_finalize
