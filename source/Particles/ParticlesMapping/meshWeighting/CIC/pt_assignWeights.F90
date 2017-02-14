!!****if* source/Particles/ParticlesMapping/meshWeighting/CIC/pt_assignWeights
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


  if (fineCoarseBdry.eqv..false.) then

     if (h(IAXIS)<0.) then
        wx(CENTER) = (1.+h(IAXIS))
        wx(LEFT_EDGE) =-h(IAXIS)
        wx(RIGHT_EDGE) = 0.
     else 
        wx(CENTER) = (1.-h(IAXIS))
        wx(LEFT_EDGE) = 0.
        wx(RIGHT_EDGE) = h(IAXIS)
     end if

     if (h(JAXIS)<0.) then
        wy(CENTER) = (1.+h(JAXIS))
        wy(LEFT_EDGE) =-h(JAXIS)
        wy(RIGHT_EDGE) = 0.
     else 
        wy(CENTER) = (1.-h(JAXIS))
        wy(LEFT_EDGE) = 0.
        wy(RIGHT_EDGE) = h(JAXIS)
     end if
     
     if (h(KAXIS)<0.) then
        wz(CENTER) = (1.+h(KAXIS))
        wz(LEFT_EDGE) =-h(KAXIS)
        wz(RIGHT_EDGE) = 0.
     else
        wz(CENTER) = (1.-h(KAXIS))
        wz(LEFT_EDGE) = 0.
        wz(RIGHT_EDGE) = h(KAXIS)
     end if
     
  else
     
     call Driver_abortFlash("pt_assignWeights: Haven't yet tested experimental fine-coarse boundary weights")

     wx(CENTER) = 0.5
     wx(LEFT_EDGE) = 0.25-0.5*h(IAXIS)
     wx(RIGHT_EDGE) = 0.25+0.5*h(IAXIS)
     wy(CENTER) = 0.5
     wy(LEFT_EDGE) = 0.25-0.5*h(JAXIS)
     wy(RIGHT_EDGE) = 0.25+0.5*h(JAXIS)
     wz(CENTER) = 0.5
     wz(LEFT_EDGE) = 0.25-0.5*h(KAXIS)
     wz(RIGHT_EDGE) = 0.25+0.5*h(KAXIS)
     
  endif
  
  return
end subroutine pt_assignWeights

