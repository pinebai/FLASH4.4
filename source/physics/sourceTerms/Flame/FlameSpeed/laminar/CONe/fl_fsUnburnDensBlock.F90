!!****if* source/physics/sourceTerms/Flame/FlameSpeed/laminar/CONe/fl_fsUnburnDensBlock
!!
!! NAME
!!
!!  fl_fsUnburnDensBlock
!!
!! SYNOPSIS
!!
!!  call fl_fsUnburnDensBlock(real,dimension(:,:,:,:),pointer  :: solndata,
!!                            real,dimension(:,:,:),intent(OUT)  :: dens,
!!                            integer,dimension(LOW:HIGH,MDIM), intent(IN)  :: calclimits)
!!
!! DESCRIPTION
!!
!!
!! Aaron Jackson, Dean Townsley, Alan Calder 2008
!!
!! ARGUMENTS
!!
!!   solndata : 
!!
!!   dens : 
!!
!!   calclimits : 
!!
!!
!!
!!***

subroutine fl_fsUnburnDensBlock(solnData, dens, calcLimits)

  use Driver_interface, ONLY: Driver_abortFlash
  use fl_fsLaminarInterface, ONLY: fl_fsConeInterp, fl_fsUnburnDens

  implicit none
#include "constants.h"
#include "Flash.h"

  real,dimension(:,:,:,:),pointer :: solnData
  real,dimension(:,:,:),intent(OUT) :: dens
  integer,dimension(LOW:HIGH,MDIM), intent(IN) :: calcLimits

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
        enddo
     enddo
  enddo

  return
end subroutine fl_fsUnburnDensBlock

!------------------------------------------------------------------------
