!!****if* source/Particles/localAPI/pt_assignWeights
!! NAME
!!
!!  pt_assignWeights
!!
!! SYNOPSIS
!!  pt_assignWeights( logical(in) :: fineCoarseBdry,
!!                    real(in)  :: h(1:MDIM),
!!                    real(out) :: wx(LEFT_EDGE:RIGHT_EDGE),
!!                    real(out) :: wy(LEFT_EDGE:RIGHT_EDGE),
!!                    real(out) :: wz(LEFT_EDGE:RIGHT_EDGE))
!!
!!
!! DESCRIPTION
!!   Computes weight factors for the fraction of charge that will be assigned
!!   to neighboring mesh cells based on the Cloud-In-Cell (CIC) scheme.
!!
!! ARGUMENTS
!!   fineCoarseBdry - Indicates if the block is a fine one at a fine-coarse boundary
!!   h              -  some sort of distance of coordinate
!!   wx             -  weights in the x direction
!!   wy             -  weights in the y direction
!!   wz             -  weights in the z direction
!!
!!  
!!***

subroutine pt_assignWeights(fineCoarseBdry,h,wx,wy,wz)

#include "constants.h"
  
  use Driver_interface, ONLY : Driver_abortFlash

  implicit none
  logical,intent(IN) :: fineCoarseBdry
  real,dimension(MDIM), intent(IN) :: h
  real,dimension(LEFT_EDGE:RIGHT_EDGE), intent(OUT) :: wx,wy,wz

  wx(:)=0.0
  wy(:)=0.0
  wz(:)=0.0

RETURN
END SUBROUTINE pt_assignWeights

