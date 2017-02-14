!!****if* source/physics/TreeRay/TreeRayMain/OpticalDepth/tr_odInit
!!
!! NAME
!!
!!  tr_odInit
!!  
!! SYNOPSIS
!!
!!  tr_odInit()
!!
!! DESCRIPTION
!!
!!
!! ARGUMENTS
!!
!!  
!!
!! PARAMETERS
!! 
!!
!! NOTES
!!   
!!
!!***

subroutine tr_odInit ()

  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
  use Driver_interface, ONLY : Driver_abortFlash
  use TreeRay_data, ONLY : tr_bhUseTreeRay, tr_nSide, tr_nPix, tr_nCd, &
    tr_meshMe
  use tr_odData, ONLY : tr_odMaxDist, tr_odCDTOIndex, tr_odICDTO, &
    tr_odICDH2, tr_odICDCO

  implicit none
#include "constants.h"
#include "Flash.h"

  integer :: istat
   

!==============================================================================

  call RuntimeParameters_get ("tr_odMaxDist", tr_odMaxDist)
  call RuntimeParameters_get ("tr_odCDTOIndex", tr_odCDTOIndex)

  ! allocate array with column densities in all directions for a single block
  tr_nCd = 1
  tr_odICDTO = 1
  if (tr_meshMe == MASTER_PE) &
  & print *, "clInit: ", tr_nCd, tr_odICDTO

#if CHEMISTRYNETWORK == 4 || CHEMISTRYNETWORK == 5
  tr_nCd = tr_nCd + 1
  tr_odICDH2 = tr_nCd
#endif
#if CHEMISTRYNETWORK == 5
  tr_nCd = tr_nCd + 1
  tr_odICDCO = tr_nCd
#endif

end subroutine tr_odInit
