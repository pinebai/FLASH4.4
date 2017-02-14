!!****if* source/physics/sourceTerms/Flame/FlameSpeed/laminar/CONe/fl_fsLaminarFlameSpeedBlock
!!
!! NAME
!!
!!  fl_fsLaminarFlameSpeedBlock
!!
!! SYNOPSIS
!!
!!  call fl_fsLaminarFlameSpeedBlock(real,dimension(:,:,:,:),pointer  :: solndata,
!!                                   real,dimension(:,:,:),intent(OUT)  :: s,
!!                                   real,dimension(:,:,:),intent(OUT)  :: dens,
!!                                   integer,dimension(LOW:HIGH,MDIM), intent(IN)  :: calclimits,
!!                                   real,OPTIONAL,dimension(:,:,:),intent(OUT)  :: ds)
!!
!! DESCRIPTION
!!
!! Aaron Jackson, Dean Townsley, Alan Calder 2008
!!
!! ARGUMENTS
!!
!!   solndata : 
!!
!!   s : 
!!
!!   dens : 
!!
!!   calclimits : 
!!
!!   ds : 
!!
!!
!!
!!***




subroutine fl_fsLaminarFlameSpeedBlock(solnData, s, dens, calcLimits, ds)

  use Driver_interface, ONLY: Driver_abortFlash
  use fl_fsLaminarInterface, ONLY: fl_fsConeInterp, fl_fsUnburnDens

  implicit none
#include "constants.h"
#include "Flash.h"

  real,dimension(:,:,:,:),pointer :: solnData
  real,dimension(:,:,:),intent(OUT) :: s
  real,dimension(:,:,:),intent(OUT) :: dens
  integer,dimension(LOW:HIGH,MDIM), intent(IN) :: calcLimits
  real,OPTIONAL,dimension(:,:,:),intent(OUT) :: ds

  real, parameter    :: ye12 = 0.5          ! ye for C12
  real, parameter    :: ye16 = 0.5          ! ye for O16
  real, parameter    :: ye22 = 10.0/22.0    ! ye for Ne22
  real :: ye                                ! ye
  real :: pres, c12i, ne22i

  integer :: i,j,k

  ! infer unburned density from pressure
  do k = calcLimits(LOW,KAXIS), calcLimits(HIGH,KAXIS)
     do j = calcLimits(LOW,JAXIS), calcLimits(HIGH,JAXIS)
        do i = calcLimits(LOW,IAXIS), calcLimits(HIGH,IAXIS)
           ! this correspondence is good for dens >~ 1.5e5 (degenerate)
           ! 0.92 is adjusted for fit
           ! constants for non-rel and relativistic degenerate e^- gas from Hansen & Kawaler
           pres = solnData(PRES_VAR,i,j,k)
#ifdef CI_MSCALAR
           c12i = solnData(CI_MSCALAR,i,j,k)
#else
           c12i = 0.5
#endif
#ifdef NEI_MSCALAR
           ne22i = solnData(NEI_MSCALAR,i,j,k)
#else
           ne22i = 0.0
#endif
           ye=c12i*ye12 + ne22i*ye22 + (1.0-c12i-ne22i)*ye16
           call fl_fsUnburnDens( pres, ye, dens(i,j,k) )
           call fl_fsConeInterp(c12i,ne22i,dens(i,j,k),s(i,j,k),ds(i,j,k))
        enddo
     enddo
  enddo

  return
end subroutine fl_fsLaminarFlameSpeedBlock

!------------------------------------------------------------------------
