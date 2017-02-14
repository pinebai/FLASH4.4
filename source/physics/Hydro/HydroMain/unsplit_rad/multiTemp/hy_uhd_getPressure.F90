!!****if* source/physics/Hydro/HydroMain/unsplit_rad/multiTemp/hy_uhd_getPressure
!!
!! NAME
!!
!!  hy_uhd_getPressure
!!
!!
!! SYNOPSIS
!!
!!  call hy_uhd_getPressure(real, intent(in)  :: eele,
!!                          real, intent(in)  :: eion,
!!                          real, intent(in)  :: erad,
!!                          real, intent(in)  :: soln(NUNK_VARS),
!!                          real, intent(out) :: pele,
!!                          real, intent(out) :: pion,
!!                          real, intent(out) :: prad) 
!!
!! DESCRIPTION
!! 
!! Calculates pele, ion and rad given energies and solution vector
!!
!! ARGUMENTS
!!
!!  eele        - electron energy 
!!  eion        - ion energy 
!!  erad        - radiation energy 
!!  soln        - solution vector     
!!  pele        - electron pressure
!!  pion        - ion pressure
!!  prad        - radiation pressure
!!
!!
!! PARAMETERS
!!
!!
!!
!!***

subroutine hy_uhd_getPressure( &
     eele, eion, erad, &
     soln, &
     pele, pion, prad)

#include "Flash.h"
#include "constants.h"
#include "Eos.h"

  use Eos_interface, ONLY: Eos

  implicit none

  ! Arguments:
  real, intent(in)  :: eele
  real, intent(in)  :: eion
  real, intent(in)  :: erad
  real, intent(in)  :: soln(NUNK_VARS)
  real, intent(out) :: pele
  real, intent(out) :: pion
  real, intent(out) :: prad
  
  ! Local variables:
  integer :: n

#ifndef FLASH_MULTISPECIES
#ifdef SUMY_MSCALAR
  integer,parameter :: sumy_map = SUMY_MSCALAR
#endif
#ifdef YE_MSCALAR
  integer,parameter :: ye_map = YE_MSCALAR
#endif
#endif

  integer :: mode, vecLen
  real    :: massfrac(NSPECIES)
  real    :: eos_arr(EOS_NUM)
  logical :: mask(EOS_VARS+1:EOS_NUM)

  ! Subroutine body:
  vecLen = 1

  ! Load inputs into EOS array. Since we are calling with
  ! MODE_DENS_EI_GATHER
  eos_arr = 0.0
  eos_arr(EOS_DENS) = soln(DENS_VAR)
  eos_arr(EOS_EINTELE) = eele
  eos_arr(EOS_EINTION) = eion
  eos_arr(EOS_EINTRAD) = erad
  eos_arr(EOS_TEMP)    = soln(TEMP_VAR) !Some Eos implementations like to find something nonzero here
  eos_arr(EOS_TEMPELE) = soln(TELE_VAR)
  eos_arr(EOS_TEMPION) = soln(TION_VAR)
  eos_arr(EOS_TEMPRAD) = soln(TRAD_VAR)
#ifndef FLASH_MULTISPECIES
#ifdef SUMY_MSCALAR
  if (soln(sumy_map).NE.0.0) then
     eos_arr(EOS_ABAR) =  1.0 /  soln(sumy_map)
  else
     eos_arr(EOS_ABAR) =  huge(soln(sumy_map))
  end if
#ifdef YE_MSCALAR
  if (soln(sumy_map).NE.0.0) then
     eos_arr(EOS_ZBAR) = soln(ye_map) /  soln(sumy_map)
  else if (soln(ye_map)==0.0) then
     eos_arr(EOS_ZBAR) = 0.0
  else
     eos_arr(EOS_ZBAR) = huge(soln(ye_map))
  end if
#endif
#endif
#endif

  ! Use the mask array to ask EOS for pressures:
  mask = .false.
  mask(EOS_PRESELE) = .true.
  mask(EOS_PRESION) = .true.
  mask(EOS_PRESRAD) = .true.

  ! Load species, if defined. Checking avoids annoying compiler
  ! warnings:
#if (NSPECIES > 0)
  do n = 1, NSPECIES
     massfrac(n) = soln(SPECIES_BEGIN-1+n)
  enddo
#endif

  call Eos(MODE_DENS_EI_COMP,vecLen,eos_arr,massfrac,mask)

  ! Recover the output:
  pele = eos_arr(EOS_PRESELE)
  pion = eos_arr(EOS_PRESION)
  prad = eos_arr(EOS_PRESRAD)

end subroutine hy_uhd_getPressure
