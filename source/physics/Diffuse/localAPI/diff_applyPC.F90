!!****if* source/physics/Diffuse/localAPI/diff_applyPC
!!
!!  NAME 
!!
!!  diff_applyPC
!!
!!  SYNOPSIS
!!
!!  call diff_applyPC(Z, R, NZ, N, LU, JA, IA, UPTR, blkLimits, blkLimitsGC)
!!
!!  DESCRIPTION 
!!
!!
!! ARGUMENTS
!!
!!
!! SIDE EFFECTS
!!
!!  
!! NOTES:
!!  
!!
!!***

subroutine diff_applyPC (Z, R, NZ, N, LU, JA, IA, UPTR, blkLimits, blkLimitsGC) 
  
  
  implicit none

#include "constants.h"

  integer, dimension(LOW:HIGH,MDIM),intent(IN) :: blkLimits, blkLimitsGC
  
  integer, intent (IN) :: NZ
  integer, intent (IN) :: N   
  real,    intent (IN) :: LU(NZ)
  integer,    intent (IN) :: JA(NZ)
  integer,    intent (IN) :: IA(N+1)
  integer,    intent (IN) :: UPTR(N)
  
  real, intent(OUT)    :: Z(blkLimitsGC(LOW,IAXIS):blkLimitsGC(HIGH,IAXIS), &
                            blkLimitsGC(LOW,JAXIS):blkLimitsGC(HIGH,JAXIS), &
                            blkLimitsGC(LOW,KAXIS):blkLimitsGC(HIGH,KAXIS))  
  real, intent(IN)     :: R(blkLimitsGC(LOW,IAXIS):blkLimitsGC(HIGH,IAXIS), &
                            blkLimitsGC(LOW,JAXIS):blkLimitsGC(HIGH,JAXIS), &
                            blkLimitsGC(LOW,KAXIS):blkLimitsGC(HIGH,KAXIS))

  
  
  Z = 0.0  
  
end subroutine diff_applyPC
