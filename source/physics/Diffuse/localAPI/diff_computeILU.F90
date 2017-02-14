!!****if* source/physics/Diffuse/localAPI/diff_computeILU
!!
!!  NAME 
!!
!!  diff_computeILU
!!
!!  SYNOPSIS
!!
!!  call diff_computeILU(AA, JA, IA, UPTR, NZ)
!!
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

!!REORDER(4): solnVec

subroutine diff_computeILU(AA, JA, IA, UPTR, NZ, N)    

  
  implicit none
  
#include "Flash.h"
#include "constants.h"  
  
  integer, intent (IN)     :: NZ
  integer, intent (IN)     :: N 
  real,    intent (INOUT)  :: AA  (NZ)
  integer,    intent (IN)  :: JA  (NZ)
  integer,    intent (IN)  :: IA  (N+1)
  integer,    intent (IN ) :: UPTR(N)

    
  
end subroutine Diff_computeILU
