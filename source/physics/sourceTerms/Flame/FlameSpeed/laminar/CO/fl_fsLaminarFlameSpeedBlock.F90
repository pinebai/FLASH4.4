!!****if* source/physics/sourceTerms/Flame/FlameSpeed/laminar/CO/fl_fsLaminarFlameSpeedBlock
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
!! Dean Townsley 2008
!! This is Shimon Asida's fit of Laminar flame speed as function of fuel density
!! From data tabulated in Timmes & Woosley 1992, ApJ, 396, 649
!! Fit only for carbon fraction Xc = 0.5 
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

  implicit none
#include "constants.h"
#include "Flash.h"

  real,dimension(:,:,:,:),pointer :: solnData
  real,dimension(:,:,:),intent(OUT) :: s
  real,dimension(:,:,:),intent(OUT) :: dens
  integer,dimension(LOW:HIGH,MDIM), intent(IN) :: calcLimits
  real,OPTIONAL,dimension(:,:,:),intent(OUT) :: ds !UNUSED in this implementation

  real, parameter :: yei = 2.0                         ! ye inverse
  real :: pres, ldens

  real, parameter    :: c1 = -43.0e0
  real, parameter    :: c2 =  4.534e0
  real, parameter    :: c3 = -0.08333e0

  integer :: i,j,k

  ! infer unburned density from pressure
  do k = calcLimits(LOW,KAXIS), calcLimits(HIGH,KAXIS)
     do j = calcLimits(LOW,JAXIS), calcLimits(HIGH,JAXIS)
        do i = calcLimits(LOW,IAXIS), calcLimits(HIGH,IAXIS)
           ! this correspondence is good for dens >~ 1.5e5 (degenerate)
           ! 0.92 is adjusted for fit
           ! constants for non-rel and relativistic degenerate e^- gas from Hansen & Kawaler
           pres = solnData(PRES_VAR,i,j,k)
           dens(i,j,k) = 0.92*yei*sqrt( (pres/1.243e15)**(6.0/4.0) + &
                                      (pres/1.004e13)**(6.0/5.0) )
           ldens = log(dens(i,j,k))
           s(i,j,k) = exp(c1 + c2*ldens + c3*(ldens**2))
        enddo
     enddo
  enddo

  return
end subroutine fl_fsLaminarFlameSpeedBlock

!------------------------------------------------------------------------
