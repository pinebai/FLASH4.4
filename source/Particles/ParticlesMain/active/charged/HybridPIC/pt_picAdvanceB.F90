!!****if* source/Particles/ParticlesMain/active/charged/HybridPIC/pt_picAdvanceB
!!
!! NAME
!!
!!  pt_picAdvanceB
!!
!! SYNOPSIS
!!
!!  call pt_picAdvanceB(real(in) :: dt)
!!
!! DESCRIPTION
!!
!!  Advance the grid magnetic field in time by dt
!!
!! ARGUMENTS
!!
!!   dt : delta t
!!
!!
!!
!!***

subroutine pt_picAdvanceB(dt)
  use Grid_interface, ONLY : Grid_addToVar
  use pt_picData, only: pt_picNsub
  use pt_picInterface, ONLY : pt_picEfield, pt_picCurl
#include "Flash.h"
#include "constants.h"  
#include "Particles.h"

  implicit none
  real, intent(in) :: dt
  real    :: h
  integer :: p, m

  m = pt_picNsub

  ! Copy B -> B0
  call Grid_addToVar(GRBX_VAR, GBX0_VAR, 1.0, .true.)
  call Grid_addToVar(GRBY_VAR, GBY0_VAR, 1.0, .true.)
  call Grid_addToVar(GRBZ_VAR, GBZ0_VAR, 1.0, .true.)

  if (mod(m,2) == 0) m = m+1   ! m must be odd
  h = dt/m 

  ! Advance B by h
  call pt_picEfield(GRBX_VAR, GRBY_VAR, GRBZ_VAR)
  call pt_picCurl(GREX_VAR, GREY_VAR, GREZ_VAR, &
       GRBX_VAR, GRBY_VAR, GRBZ_VAR, -h, .false.)

  do p = 1, (m-1)/2
     ! Advance B0 by 2h, using B to compute E
     call pt_picEfield(GRBX_VAR, GRBY_VAR, GRBZ_VAR)
     call pt_picCurl(GREX_VAR, GREY_VAR, GREZ_VAR, &
          GBX0_VAR, GBY0_VAR, GBZ0_VAR, -2*h, .false.)
     ! Advance B by 2h, using B0 to compute E
     call pt_picEfield(GBX0_VAR, GBY0_VAR, GBZ0_VAR)
     call pt_picCurl(GREX_VAR, GREY_VAR, GREZ_VAR, &
          GRBX_VAR, GRBY_VAR, GRBZ_VAR, -2*h, .false.)
  end do

  ! Advance B0 by h
  call pt_picEfield(GRBX_VAR, GRBY_VAR, GRBZ_VAR)
  call pt_picCurl(GREX_VAR, GREY_VAR, GREZ_VAR, &
       GBX0_VAR, GBY0_VAR, GBZ0_VAR, -h, .false.)

  ! B and B0 at time p*h.  Average B and B0, store in B
  call Grid_addToVar(GRBX_VAR, GRBX_VAR, -0.5, .false.) ! 0.5*B
  call Grid_addToVar(GRBY_VAR, GRBY_VAR, -0.5, .false.)
  call Grid_addToVar(GRBZ_VAR, GRBZ_VAR, -0.5, .false.)
  
  call Grid_addToVar(GBX0_VAR, GRBX_VAR, 0.5, .false.) ! + 0.5*B0
  call Grid_addToVar(GBY0_VAR, GRBY_VAR, 0.5, .false.)
  call Grid_addToVar(GBZ0_VAR, GRBZ_VAR, 0.5, .false.)

end subroutine pt_picAdvanceB
