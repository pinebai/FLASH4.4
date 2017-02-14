!!***if* source/Particles/ParticlesMapping/meshWeighting/pt_mapOneParticle
!!
!! NAME
!!
!!  pt_mapOneParticle
!!
!!
!! SYNOPSIS
!!
!!  pt_mapOneParticle(integer, dimension(LOW:HIGH,MDIM), intent(IN) :: blkLimitsGC, &
!!                    integer, dimension(MDIM), intent(IN) :: intPos, &
!!                    real, intent(IN)    :: pt_attValue, &
!!                    logical, intent(IN) :: fcBdry, &
!!                    real, dimension(MDIM),intent(IN) :: hfunc
!!                    real, dimension(blkLimitsGC(LOW,IAXIS):blkLimitsGC(HIGH,IAXIS), &
!!                                    blkLimitsGC(LOW,JAXIS):blkLimitsGC(HIGH,JAXIS), &
!!                                    blkLimitsGC(LOW,KAXIS):blkLimitsGC(HIGH,KAXIS)), &
!!                             intent(INOUT) :: buff
!!
!! DESCRIPTION
!!
!! Routine which smears a single particle across the cells in the temporary 
!! buffer (buff).  The particle is smeared according to the mapping scheme 
!! specified in the subroutine pt_assignWeights.
!! 
!! ARGUMENTS
!!           blkLimitsGC:  Size of a block including guard cells.
!!           intPos:  Coordinates of the cell that holds the particle.
!!           pt_attValue:  Particle attribute that will be mapped to the mesh.
!!           fcBdry:  Specifies whether the particle smears across a fine-coarse boundary.
!!           hfunc:   Displacement of particle from cell center.  Required in pt_assignWeights.
!!           buff:  Temporary buffer in which to accumulate particle mapping.
!!
!!***

subroutine pt_mapOneParticle (blkLimitsGC,intPos,pt_attValue,fcBdry, &
     hfunc, buff)
#include "constants.h"
#include "Flash.h"

  implicit none

  integer, dimension(LOW:HIGH,MDIM), intent(IN) :: blkLimitsGC
  integer, dimension(MDIM), intent(IN) :: intPos
  real, intent(IN)    :: pt_attValue
  logical, intent(IN) :: fcBdry
  real, dimension(MDIM),intent(IN) :: hfunc
  real, dimension(blkLimitsGC(LOW,IAXIS):blkLimitsGC(HIGH,IAXIS),&
                  blkLimitsGC(LOW,JAXIS):blkLimitsGC(HIGH,JAXIS),&
                  blkLimitsGC(LOW,KAXIS):blkLimitsGC(HIGH,KAXIS)),&
                  INTENT(inout) :: buff

  real,dimension(-1:1):: wx,wy,wz
  real,dimension(-1:1,-1:1,-1:1) :: wt
  integer :: ip, jp, kp, i, j, k
  integer :: istart,iend,jstart,jend,kstart,kend

  ip = intPos(IAXIS)
  jp = intPos(JAXIS)
  kp = intPos(KAXIS)

! Pass .false. temporarily to use the traditional CIC scheme.
  call pt_assignWeights(.false.,hfunc,wx,wy,wz)

 
  istart=-1; iend=1
  if(NDIM>1) then
     jstart=-1; jend=1
  else
     jstart=0; jend=0
  end if
  if(NDIM>2) then
     kstart=-1; kend=1
  else
     kstart=0; kend=0
  end if
  

  do k = kstart,kend
     do j = jstart,jend
        do i = istart,iend
           wt(i,j,k)= wx(i)*wy(j)*wz(k)
        end do
     end do
  end do


  do k = kstart,kend
     do j = jstart,jend
        do i = istart,iend
           buff(ip+i,jp+j,kp+k) = buff(ip+i,jp+j,kp+k) + wt(i,j,k) * pt_attValue
        end do
     end do
  end do
  
  return
  
end subroutine pt_mapOneParticle

