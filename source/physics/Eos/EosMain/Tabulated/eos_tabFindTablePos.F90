!!****if* source/physics/Eos/EosMain/Tabulated/eos_tabFindTablePos
!!
!! NAME
!!
!!  eos_tabFindTablePos
!!
!! SYNOPSIS
!!
!!  call eos_tabFindTablePos(integer (in)  :: species,
!!                                    real    (in)  :: speciesTemperature,
!!                                    real    (in)  :: speciesDensity,
!!                                    integer (in)  :: maxComp,
!!                                    integer (in)  :: selTT)
!!
!! DESCRIPTION
!!
!!  Given desired Temperature and Density values and a ((temp(i), dens(k),
!!  i=1,ntemp, k=1,ndens) description of a (temperature, density) table grid,
!!  find and return the index pair (i,k) such that
!!        temp(i) < Temperature <=  temp(j)
!!  and
!!        dens(k) < Density     <=  dens(l)  .
!!
!! ARGUMENTS
!!
!!   species            : The species index
!!   speciesTemperature : The species temperature
!!   speciesDensity     : The species density
!!   maxComp :          : The highest component (ions, electrons, or matter) for which
!!                        data are to be returned
!!   selTT              : Select ZF, EN, or HC table type
!!
!!***

#include "Flash.h"
#include "Eos.h"

subroutine eos_tabFindTablePos (species,            &
                                        speciesTemperature, &
                                        speciesDensity,     &
                                        selTT,              &
                                        i,j,k,l,&
                                        tau, delta, &
                                        iPrev, kPrev, &
                                        lowerBoundaryDens, &
                                        lowerBoundaryTemp,withinBoundaryTemp,withinBoundaryDens, &
                                        T1,T2,D1,D2)

  use Driver_interface,  ONLY : Driver_abortFlash
  use eos_tabData,      ONLY : EOS_TAB_NCOMP,              &
                               EOS_TAB_NALLTAB,            &
                               EOS_TABVT_ZF,               &
                               EOS_TABVT_EN,               &
                               EOS_TABVT_HC,               &
                               EOS_TABVT_PR,               &
                               EOS_TABVT_ENTR,             &
                                eos_useLogTables,          &
                                eosT_varTableGroupPT,      &
                                eos_allTab


  implicit none

  integer, intent (in)  :: species
  real,    intent (in)  :: speciesTemperature
  real,    intent (in)  :: speciesDensity
  integer, intent (in)  :: selTT

  type(eosT_varTableGroupPT),pointer :: thisTypeTable
  real,    pointer :: thisTypeDensity (:)
  real,    pointer :: thisTypeTemperature (:)

  logical,intent(OUT) :: lowerBoundaryDens
  logical :: upperBoundaryDens
  logical,intent(OUT) :: withinBoundaryDens
  logical,intent(OUT) :: lowerBoundaryTemp
  logical :: upperBoundaryTemp
  logical,intent(OUT) :: withinBoundaryTemp
  logical :: found

  integer :: d,t
  integer,intent(OUT) :: i,j,k,l
  integer,intent(INOUT) :: iPrev,kPrev
  integer :: varType
  integer :: nstepsDensity
  integer :: nstepsTemperature
  integer :: g

  real,intent(OUT) :: D1,D2,T1,T2
  real :: f1,f2,f3,f4
  real :: o1,o2,o3,o4
  real :: speciesDensityTT
  real :: speciesTemperatureTT
  real,intent(out) :: tau,delta

!
!
!   ...Set the handle to the ZF, EN, or HC tables and associated data arrays.
!
!
  select case (selTT)
  case(EOS_TABULAR_Z)
     varType = EOS_TABVT_ZF
  case(EOS_TABULAR_E)
     varType = EOS_TABVT_EN
  case(EOS_TABULAR_C)
     varType = EOS_TABVT_HC
  case(EOS_TABULAR_P)
     varType = EOS_TABVT_PR
  case(EOS_TABULAR_S)
     varType = EOS_TABVT_ENTR
  case default
     varType = 0
  end select

  if (varType == 0) then
      call Driver_abortFlash ('[eos_tabFindTablePos] ERROR: no handle to tables')
  end if

  thisTypeTable => eos_allTab(species)%tg(varType)
  thisTypeDensity => thisTypeTable%td%Densities
  thisTypeTemperature => thisTypeTable%td%Temperatures
!
!
!   ...Get the current temperature and ion number density of the species and find:
!
!        1) the temperature/density quadrant (T1,T2,D1,D2) boundary containing
!           containing the species's temperature and density (x)
!
!        2) the table index quadrant (i,j,k,l) of the boundary, as indicated
!           in the figure below.
!
!        3) the four tabulated values (o1,o2,o3,o4)
!
!
!                    o3------------o4   T2 (j)
!                     |            |
!                     |            |
!                     |        x   |
!                     |            |
!                     |            |
!                    o1 -----------o2   T1 (i)
!
!                  D1 (k)         D2 (l)
!
!
!      In case the temperature and/or density of the species lay outside
!      the tabulated boundaries, take the corresponding boundary values
!      of the tables. The criteria as to when the species's temperature and
!      density belong within the boundary are:
!
!                         T1 =< speciesTemperature =< T2
!                         D1 =<   speciesDensity   =< D2
!
!
  if (eos_useLogTables) then
      speciesTemperatureTT = log10 (speciesTemperature)
      speciesDensityTT     = log10 (speciesDensity)
  else
      speciesTemperatureTT = speciesTemperature
      speciesDensityTT     = speciesDensity
  end if

  nstepsTemperature = thisTypeTable%td%ntemp

  lowerBoundaryTemp = speciesTemperatureTT < thisTypeTemperature (1                )
  if (lowerBoundaryTemp) then
     upperBoundaryTemp = .FALSE.
  else
     upperBoundaryTemp = speciesTemperatureTT > thisTypeTemperature (nstepsTemperature) .OR. nstepsTemperature == 1
  end if
  withinBoundaryTemp = (.not.lowerBoundaryTemp) .and. (.not.upperBoundaryTemp)

  if (withinBoundaryTemp) then
     if (iPrev > 0 .AND. iPrev < nstepsTemperature) then
        i = iPrev
        j = i + 1; 
     
        T1 = thisTypeTemperature (i)
        T2 = thisTypeTemperature (j)

        found = ((T1-speciesTemperatureTT)*(T2-speciesTemperatureTT) .LE. 0.0)
     else
        found = .FALSE.
     end if
     if (.NOT.found) then
        do t = 2,nstepsTemperature
           if (thisTypeTemperature (t) >= speciesTemperatureTT) then
              i  = t - 1
              j  = t
              T1 = thisTypeTemperature (i)
              T2 = thisTypeTemperature (j)
#ifndef FLASH_OPENMP              
              iPrev = i
#endif
              exit
           end if
        end do
     end if
  else
      if (lowerBoundaryTemp) then
          i  = 1
          j  = 1
          T1 = speciesTemperatureTT
          T2 = thisTypeTemperature (1)
      end if

      if (upperBoundaryTemp) then
          i  = nstepsTemperature
          j  = nstepsTemperature
          T1 = thisTypeTemperature (nstepsTemperature)
          T2 = speciesTemperatureTT
      end if
  end if

  nstepsDensity = thisTypeTable%td%ndens

  ! Check to see if density is off of the table:
  lowerBoundaryDens = speciesDensityTT < thisTypeDensity (1            )
  if (lowerBoundaryDens) then
     upperBoundaryDens = .FALSE.
  else
     upperBoundaryDens = speciesDensityTT > thisTypeDensity (nstepsDensity) .OR. nstepsDensity == 1
  end if
  withinBoundaryDens = (.not.lowerBoundaryDens) .and. (.not.upperBoundaryDens)

  ! if (withinBoundaryDens == .true.) -> Density is within table boundaries
  if (withinBoundaryDens) then
     if (kPrev > 0 .AND. kPrev < nstepsDensity) then
        k = kPrev
        l = k + 1

        D1 = thisTypeDensity (k)
        D2 = thisTypeDensity (l)

        found = ((D1-speciesDensityTT)*(D2-speciesDensityTT) .LE. 0.0)
     else
        found = .FALSE.
     end if
     if (.NOT.found) then
        do d = 2,nstepsDensity
           if (thisTypeDensity (d) >= speciesDensityTT) then
              k  = d - 1
              l  = d
              D1 = thisTypeDensity (k)
              D2 = thisTypeDensity (l)
#ifndef FLASH_OPENMP
              kPrev = k
#endif              
              exit
           end if
        end do
     end if
  else
      if (lowerBoundaryDens) then
          k  = 1
          l  = 1
          D1 = speciesDensityTT
          D2 = thisTypeDensity (1)
      end if

      if (upperBoundaryDens) then
          k  = nstepsDensity
          l  = nstepsDensity
          D1 = thisTypeDensity (nstepsDensity)
          D2 = speciesDensityTT
      end if
  end if

  tau   = (speciesTemperatureTT - T1) / (T2 - T1)
  delta = (speciesDensityTT     - D1) / (D2 - D1)


  return
end subroutine eos_tabFindTablePos
