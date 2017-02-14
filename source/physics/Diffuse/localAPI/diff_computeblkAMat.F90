!!****if* source/physics/Diffuse/localAPI/diff_computeblkAMat
!!
!!  NAME 
!!
!!  diff_computeblkAMat
!!
!!  SYNOPSIS
!!
!!  call diff_computeblkAMat(AA, JA, IA, UPTR, NZ, N, blockID, iFactorA, iFactorB, iFactorC)!!
!!
!!  DESCRIPTION 
!!
!!
!!
!!
!! ARGUMENTS
!!
!!
!! SIDE EFFECTS
!!
!!  
!! NOTES
!!
!!
!!***

!!REORDER(4): solnVec

subroutine diff_computeblkAMat (blockID, NZ, N, dt, theta, AA, JA, IA, UPTR, iFactorA, iFactorB, iFactorC)
  
  
  implicit none
  
#include "Flash.h"
#include "constants.h"  
  
  integer, intent (IN)  :: blockID
  integer, intent (IN)  :: NZ
  integer, intent (IN)  :: N 
  real,    intent (IN)  :: dt
  real,    intent (IN)  :: theta  
  real,    intent (OUT) :: AA  (NZ)
  integer, intent (OUT) :: JA  (NZ)
  integer, intent (OUT) :: IA  (N+1)
  integer, intent (OUT) :: UPTR(N)  
  integer, intent(IN)   :: iFactorA
  integer, intent(IN)   :: iFactorB
  integer, OPTIONAL, intent(IN) :: iFactorC  
  
  AA(:) = 0.0
  JA(:) = 0
  IA(:) = 0
  UPTR(:) = 0
  
end subroutine diff_computeblkAMat
