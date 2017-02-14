!!****if* source/physics/Diffuse/localAPI/diff_computeAX
!!
!!  NAME 
!!
!!  diff_computeAX
!!
!!  SYNOPSIS
!!
!!  call diff_computeAX (integer(IN)           :: blkLimits,
!!                  integer(IN)           :: blkLimitsGC,
!!                  real(IN)              :: X    ,
!!                  real(OUT)             :: AX ,
!!                  real(IN)              :: iFactorA    ,
!!                  real(IN)   , OPTIONAL :: iFactorC ,
!!                  real(IN)   , OPTIONAL :: iFactorD ,
!!                  real(IN)              :: dt,
!!                  real(IN)              :: theta)
!!
!!
!!  DESCRIPTION 
!!      This routine computes AX, Jacobian free, as we do not have to store AX.
!!      The elements of A are built from FD of general diffusion/conduction equation      
!!
!!      A * dV/dt = d/dx(B*dV/dx) + d/dy(B*dV/dy) + d/dx(B*dV/dz) + C*V + D
!!
!!
!! ARGUMENTS
!!      
!!  blkLimits - integer index limits for the interior of the block
!!  blkLimitsGC -integer index limits for the whole block
!!  X - input vector
!!  AX - computed value to return
!!  iFactorA -    
!!  iFactorB -    
!!  iFactorC -    
!!  del -
!!  dt - 
!!  theta - 
!!
!!      
!!
!! SIDE EFFECTS
!!      
!!  
!! NOTES:
!!  
!!
!!***

subroutine diff_computeAX (blockID, blkLimits, blkLimitsGC, iVar, iFactorA, dt, theta, AX,iFactorC, iFactorD)
  implicit none 
  
#include "constants.h"
  
  !!------------------------------------------------------------------------------------------
  integer, intent(IN) :: blockID
  integer, dimension(2,MDIM),intent(IN) :: blkLimits
  integer, dimension(2,MDIM),intent(IN) :: blkLimitsGC
  integer,intent(IN) :: iVar
  integer,intent(IN) :: iFactorA
  real,intent(IN):: dt
  real,intent(IN):: theta   
  real,intent(OUT) :: AX (blkLimitsGC(LOW,IAXIS):blkLimitsGC(HIGH,IAXIS), &
       blkLimitsGC(LOW,JAXIS):blkLimitsGC(HIGH,JAXIS), &
       blkLimitsGC(LOW,KAXIS):blkLimitsGC(HIGH,KAXIS))
  integer,intent(IN),OPTIONAL :: iFactorC
  integer,intent(IN),OPTIONAL :: iFactorD

  AX(:,:,:) = 0.0

  !!-------------------------------------------------------------------------------------------
  
end subroutine diff_computeAX
