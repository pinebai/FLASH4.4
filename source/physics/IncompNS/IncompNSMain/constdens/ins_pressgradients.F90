!!****if* source/physics/IncompNS/IncompNSMain/constdens/ins_pressgradients
!!
!! NAME
!!
!!  ins_pressgradients
!!
!! SYNOPSIS
!!
!!  call ins_pressgradients(real(in) :: time,
!!                          real(in) :: dt)
!!
!! DESCRIPTION
!!
!!
!! ARGUMENTS
!!
!!   time : 
!!
!!   dt : 
!!
!! AUTOGENROBODOC
!!
!!
!!***


#include "Flash.h"
#include "constants.h"
#include "IncompNS.h"


subroutine ins_pressgradients(time,dt)

  use IncompNS_data, only : ins_meshMe,ins_dpdx,ins_dpdy,ins_dpdz, &
                            ins_constmass,ins_WB,ins_WBold,ins_WBREF
  use ins_interface, only : ins_getBulkVelocity

  implicit none
  !! ---- Argument List ----------------------------------
  real, intent(in) :: time, dt  
  !! -----------------------------------------------------

  ! Local variables:
#ifdef EXPONENTIAL_RAMP
  real, parameter :: tau=0.025 
#endif 
  real :: rampfact, WBREF

  ! No need to assign gradsP fixed on runtime parameters file, as they
  ! have been loaded in IncompNS_init.

  ! Run at constant mass in Z?
  if (ins_constmass) then

    call ins_getBulkVelocity(ins_WB,KAXIS)

    ! Estimate DPDZ:
    ! Ramp factor:
    rampfact = 1.
#ifdef EXPONENTIAL_WBREF_RAMP
    !if (ins_WB .lt. ins_WBREF) 
    rampfact = (rampfact-exp(-time/tau))
#endif
    WBREF = rampfact*ins_WBREF
    ins_dpdz = ins_dpdz + ( 1.1*(ins_WB - WBREF)/dt - 1.*(ins_WBold - WBREF)/dt )

    if (ins_meshMe .eq. MASTER_PE) then
      print*, 'Constmass, WB=',ins_WB,', WBold=',ins_WBold,', dpdz=',ins_dpdz
      print*, 'Reference  WB=',ins_WBREF,', WBrmp=',WBREF,', dti=',dt
    endif

    ins_WBold = ins_WB

  else

    !if (ins_meshMe .eq. MASTER_PE) print*, 'dpdx=',ins_dpdx,', dpdy=',ins_dpdy,', dpdz=',ins_dpdz

  endif

  return

end subroutine 
