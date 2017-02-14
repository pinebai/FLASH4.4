!!****if* source/physics/sourceTerms/Flame/FlameSpeed/BuoyancyCompensation/CONe/fl_fsAtwood
!!
!! NAME
!!
!!  fl_fsAtwood
!!
!! SYNOPSIS
!!
!!  call fl_fsAtwood(real, dimension(:,:,:,:), pointer  :: solndata,
!!                   real, dimension(:,:,:)(out) :: atwood,
!!                   real, dimension(:,:,:)(in) :: dens_u,
!!                   integer, dimension(LOW:HIGH,MDIM)(in) :: calclimits)
!!
!! DESCRIPTION
!!
!! fill the atwood array with an estimate of the atwood number
!! for fuel starting at density dens_u
!! the actual atwood number as a function of density has been
!! calculated and tabulated at startup
!! for densities outside of the tabulated range we extrapolate
!!
!! Aaron Jackson, Dean Townsley, Alan Calder 2008
!!
!! ARGUMENTS
!!
!!   solndata : solution data
!!
!!   atwood : atwood array
!!
!!   dens_u : density u
!!
!!   calclimits : calculation limits
!!
!!
!!***


subroutine fl_fsAtwood(solnData, atwood, dens_u, calcLimits)

  use fl_fsAtwoodData, ONLY : fl_fsAtwoodTabMinLdens, fl_fsAtwoodTabMaxLDens, &
                              fl_fsAtwoodTabDLdens, fl_fsAtwoodTabLdensNum,   &
                              fl_fsAtwoodTabMinXc12, fl_fsAtwoodTabMaxXc12,   &
                              fl_fsAtwoodTabDXc12, fl_fsAtwoodTabXc12Num,     &
                              fl_fsAtwoodTabA
  use Driver_interface, ONLY : Driver_abortFlash
#include "constants.h"
#include "Flash.h"
  implicit none

  real, dimension(:,:,:,:), pointer :: solnData
  real, dimension(:,:,:), intent(out) :: atwood
  real, dimension(:,:,:), intent(in)  :: dens_u
  integer, dimension(LOW:HIGH,MDIM), intent(in)    :: calcLimits

  integer :: i,j,k
  integer :: tab_i, tab_j
  real    :: x1, x2, ldens, c12
  real    :: atwood11,atwood12,atwood21,atwood22
  real    :: c1,c2  !! interpolation constants

  do k = calcLimits(LOW,KAXIS), calcLimits(HIGH,KAXIS)
     do j = calcLimits(LOW,JAXIS), calcLimits(HIGH,JAXIS)
        do i = calcLimits(LOW,IAXIS), calcLimits(HIGH,IAXIS)
           ldens = log(dens_u(i,j,k))
           c12   = solnData(CI_MSCALAR,i,j,k)
           ! "coordinate" in table grid
           tab_i = int((ldens-fl_fsAtwoodTabMinLdens)/fl_fsAtwoodTabDLdens+1)
           tab_j = int((c12-fl_fsAtwoodTabMinXc12)/fl_fsAtwoodTabDXc12+1)
           ! extrapolate if outside table
           if (tab_i < 1) tab_i = 1
           if (tab_j < 1) tab_j = 1
           ! note table has Num+1 entries (Num counts intervals)
           if (tab_i > fl_fsAtwoodTabLdensNum) tab_i = fl_fsAtwoodTabLdensNum
           if (tab_j > fl_fsAtwoodTabXc12Num) tab_j = fl_fsAtwoodTabXc12Num
           x1 = fl_fsAtwoodTabMinLdens + (tab_i-1)*fl_fsAtwoodTabDLdens
           x2 = fl_fsAtwoodTabMinXc12 + (tab_j-1)*fl_fsAtwoodTabDXc12
           c1 = (ldens - x1)/fl_fsAtwoodTabDLdens
           c2 = (c12 - x2)/fl_fsAtwoodTabDXc12
           atwood11=fl_fsAtwoodTabA(tab_i,tab_j)
           atwood21=fl_fsAtwoodTabA(tab_i+1,tab_j)
           atwood12=fl_fsAtwoodTabA(tab_i,tab_j+1)
           atwood22=fl_fsAtwoodTabA(tab_i+1,tab_j+1)
           atwood(i,j,k) = (1.0-c1)*(1.0-c2)*atwood11 + c1*(1.0-c2)*atwood21 &
                         + c1*c2*atwood22 + (1.0-c1)*c2*atwood12
        enddo
     enddo
  enddo

end subroutine
