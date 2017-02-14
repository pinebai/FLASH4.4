!!****if* source/physics/Eos/EosMain/Tabulated/eos_tabUpdateOutsideCount
!!
!! NAME
!!
!!  eos_tabUpdateOutsideCount
!!
!! SYNOPSIS
!!
!!  call eos_tabUpdateOutsideCount(integer(in) :: species,
!!                                 logical(in) :: upperboundarytemp,
!!                                 logical(in) :: tempislog,
!!                                 real(in) :: temp,
!!                                 logical(in) :: upperboundarydens,
!!                                 logical(in) :: densislog,
!!                                 real(in) :: dens)
!!
!! DESCRIPTION
!!
!!
!! ARGUMENTS
!!
!!   species : species
!!
!!   upperboundarytemp : check if upper boundary for temperature
!!
!!   tempislog : check if temperature is log
!!
!!   temp : temperature 
!!
!!   upperboundarydens : check if upper boundary for density
!!
!!   densislog : check if density is log
!!
!!   dens : density
!!
!!
!!
!!***

#include "Eos.h"

subroutine eos_tabUpdateOutsideCount(species, &
                                     upperBoundaryTemp, tempIsLog, &
                                     temp, &
                                     upperBoundaryDens, densIsLog, &
                                     dens)

  use Eos_data,                   ONLY : eos_logLevel
  use eos_tabData,                ONLY : EOS_TAB_NALLTAB,           &
                                         eos_allTab, &
                                         eos_tabAllDiag
  implicit none

  integer,intent(in) :: species
  logical,intent(in) :: upperBoundaryTemp, upperBoundaryDens
  logical,intent(in) :: tempIsLog, densIsLog
  real,   intent(in) :: temp, dens

  real :: physTemp, physDens


  if (eos_logLevel .LT. EOS_LOGLEVEL_WARN_ANY) then
     return                     !RETURN IMMEDIATELY
  end if
  if (.NOT. (upperBoundaryTemp .OR. upperBoundaryDens)) then
     return                     !RETURN IMMEDIATELY
  end if
  
  if (tempIsLog) then
     physTemp = 10.0**temp
  else
     physTemp = temp
  end if
  if (densIsLog) then
     physDens = 10.0**dens
  else
     physDens = dens
  end if

  if (upperBoundaryTemp) then
     if (eos_tabAllDiag(species)%highTempCount == 0) then
        eos_tabAllDiag(species)%firstHighTempEvent%temp = physTemp
        eos_tabAllDiag(species)%firstHighTempEvent%dens = physDens
     end if
     eos_tabAllDiag(species)%highTempCount = &
          eos_tabAllDiag(species)%highTempCount + 1
     if (physTemp > eos_tabAllDiag(species)%highestTemp) &
          eos_tabAllDiag(species)%highestTemp = physTemp
  end if

  if (upperBoundaryDens) then
     if (eos_tabAllDiag(species)%highDensCount == 0) then
        eos_tabAllDiag(species)%firstHighDensEvent%temp = physTemp
        eos_tabAllDiag(species)%firstHighDensEvent%dens = physDens
     end if
     eos_tabAllDiag(species)%highDensCount = &
          eos_tabAllDiag(species)%highDensCount + 1
     if (physDens > eos_tabAllDiag(species)%highestDens) &
          eos_tabAllDiag(species)%highestDens = physDens
  end if
  

end subroutine eos_tabUpdateOutsideCount
