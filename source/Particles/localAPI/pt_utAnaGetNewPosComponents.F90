!!****if* source/Particles/localAPI/pt_utAnaGetNewPosComponents
!!
!! NAME
!!
!!  pt_utAnaGetNewPosComponents
!!
!! SYNOPSIS
!!
!!  pt_utAnaGetNewPosComponents(real(inout)  :: particles(:,:),
!!                              integer(in) :: maxParticlesPerProc,
!!                              integer(in) :: xOutPart,
!!                              integer(in) :: yOutPart,
!!                              integer(in) :: zOutPart,
!!                              real(in)    :: t
!!                              integer(in) :: k)
!!
!! DESCRIPTION
!!
!!  Get new position components for the analytical solution
!!
!! ARGUMENTS
!!
!!  particles : the data structure containing the particles
!!              It is two dimensional real array, the first dimension
!!              represents properties associated with the data 
!!              structure, and 
!!              second dimension is index to individual elements in 
!!              the datastructure.
!! maxParticlesPerProc : the maximum count of particles allowed on 
!!                       any processor
!! xOutPart    : particle property with x component 
!! yOutPart    : particle property with y component 
!! zOutPart    : particle property with z component 
!! t        : time
!! k        : identity of a specific particle
!!
!!
!!***
subroutine pt_utAnaGetNewPosComponents(particles,maxParticlesPerProc,xOutPart,yOutPart,zOutPart,t,k)

  implicit none

#include "Flash.h"

  integer, INTENT(in) :: maxParticlesPerProc
  real, INTENT(inout),dimension(NPART_PROPS,maxParticlesPerProc) :: particles
  integer, INTENT(in) :: xOutPart, yOutPart, zOutPart
  real, INTENT(in)           :: t
  integer, INTENT(in)           :: k


end subroutine pt_utAnaGetNewPosComponents
