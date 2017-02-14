!!****if* source/physics/Eos/EosMain/Tabulated/eos_tabGetSpeciesTabulatedData
!!
!! NAME
!!
!!  eos_tabGetSpeciesTabulatedData
!!
!! SYNOPSIS
!!
!!  call eos_tabGetSpeciesTabulatedData (integer (in) :: species,
!!                                        real    (in) :: speciesTemperature,
!!                                        real    (in) :: speciesDensity,
!!                                        integer (in) :: maxComp,
!!                                        integer (in) :: needZFTable,
!!                                        integer (in) :: needENDerivs,
!!                                        integer (in) :: needHCDerivs)
!!
!! DESCRIPTION
!!
!!  Computes absorption, emission and transport opacities from tabulated
!!  values for a particular (temperature, density, energy group, species)
!!  quadruple.
!!
!! ARGUMENTS
!!
!!   species            : The species index
!!   speciesTemperature : The species temperature
!!   speciesDensity     : The species density
!!   needZFDerivs        : Knob to activate data extraction from the average ionization tables.
!!                         Number of highest derivative needed, may be 0, negative for none.
!!   needENDerivs        : Knob to activate data extraction from the internal energy tables
!!                         Number of highest derivative needed, may be 0, negative for none.
!!   needHCDerivs        : Knob to activate data extraction from the specific heat tables
!!                         Number of highest derivative needed, may be 0, negative for none.
!!
!! NOTES
!!
!!  See definitions of EOS_TABINT_DERIV_* in eos_tabData and EOS_TAB_NDERIVS for
!!  the numbering of derivatives.
!!***
subroutine eos_tabGetSpeciesTabulatedData (species,            &
                                            speciesTemperature, &
                                            speciesDensity,     &
                                            wantedDerivs, &
                                            outData              )

  use eos_tabInterface, ONLY: eos_tabGetSpeciesAnyTableData, eos_tabFindTablePos, &
                              eos_tabUpdateOutsideCount
  use eos_tabData, ONLY: eos_tabIonizationKind,     &
                          eos_tabIntEnergyKind,       &
                          eos_tabHeatCpKind,      &
                          EOS_TAB_NCOMP,EOS_TAB_NALLTAB, &
                          EOS_TAB_NDERIVS, &
                          EOS_TAB_FOR_ION, EOS_TAB_FOR_ELE, EOS_TAB_FOR_MAT, &
                          EOS_TABVT_ZF, EOS_TABVT_EN, EOS_TABVT_PR, EOS_TABVT_HC, &
                          EOS_TABVT_ENTR, &
                          EOS_TABINT_DERIV_0, EOS_TAB_NDERIVS,      &
                          eos_allTab
  implicit none
  
#include "Eos.h"
  
  integer, intent (in) :: wantedDerivs(EOS_TAB_NCOMP,EOS_TAB_NALLTAB)
  integer, intent (in) :: species
  real,    intent (in) :: speciesTemperature
  real,    intent (in) :: speciesDensity
  real,intent(OUT) :: outData(0:EOS_TAB_NDERIVS,1:EOS_TABVT_ENTR,EOS_TAB_NCOMP)

  logical :: lowerBoundaryDens
  logical :: withinBoundaryDens
  logical :: lowerBoundaryTemp
  logical :: withinBoundaryTemp
  logical :: isLog
  integer :: i,j,k,l
  integer :: varType
  integer,save :: iSave=1,kSave=1
  real :: D1,D2,T1,T2
  real :: tau,delta


  iSave = -1; kSave = -1
  call eos_tabFindTablePos (species,            &
                                        speciesTemperature, &
                                        speciesDensity,     &
                                        EOS_TABVT_ZF,       &
                                        i,j,k,l,&
                                        tau, delta, &
                                        iSave, kSave, &
                                        lowerBoundaryDens, &
                                        lowerBoundaryTemp,withinBoundaryTemp,withinBoundaryDens, &
                                        T1,T2,D1,D2)

  varType = EOS_TABVT_ZF
!!$  thisTypeTable => eos_allTab(species)%tg(varType)
  isLog = eos_allTab(species)%tg(varType)%td%isLog
  call eos_tabUpdateOutsideCount(species, &
                                 .NOT.(lowerBoundaryTemp .OR. withinBoundaryTemp), &
                                 isLog, T2, &
                                 .NOT.(lowerBoundaryDens .OR. withinBoundaryDens), &
                                 isLog, D2)

!
!
!   ...Extract only the necessary opacities from the tables.
!
!
  if (ANY(wantedDerivs(:,EOS_TABVT_ZF) .GE. 0)) then

      call eos_tabGetSpeciesAnyTableData (species,            &
                                        speciesTemperature, &
                                        speciesDensity,     &
                                        (wantedDerivs .GE. 0),&
                                        EOS_TABULAR_Z,      &
                                        maxval(wantedDerivs(:,EOS_TABVT_ZF)),       &
                                        i,j,k,l,&
                                        tau, delta, &
                                        lowerBoundaryDens, &
                                        lowerBoundaryTemp,withinBoundaryTemp,withinBoundaryDens, &
                                        outData(:,EOS_TABVT_ZF,:), &
                                        T1,T2,D1,D2)
  end if

  if (ANY(wantedDerivs(:,EOS_TABVT_EN) .GE. 0)) then

      call eos_tabGetSpeciesAnyTableData (species,            &
                                        speciesTemperature, &
                                        speciesDensity,     &
                                        (wantedDerivs .GE. 0),&
                                        EOS_TABULAR_E,      &
                                        maxval(wantedDerivs(:,EOS_TABVT_EN)),       &
                                        i,j,k,l,&
                                        tau, delta, &
                                        lowerBoundaryDens, &
                                        lowerBoundaryTemp,withinBoundaryTemp,withinBoundaryDens, &
                                        outData(:,EOS_TABVT_EN,:), &
                                        T1,T2,D1,D2)
  end if

  if (ANY(wantedDerivs(:,EOS_TABVT_HC) .GE. 0)) then

      call eos_tabGetSpeciesAnyTableData (species,            &
                                        speciesTemperature, &
                                        speciesDensity,     &
                                        (wantedDerivs .GE. 0),&
                                        EOS_TABULAR_C,      &
                                        maxval(wantedDerivs(:,EOS_TABVT_HC)),       &
                                        i,j,k,l,&
                                        tau, delta, &
                                        lowerBoundaryDens, &
                                        lowerBoundaryTemp,withinBoundaryTemp,withinBoundaryDens, &
                                        outData(:,EOS_TABULAR_C,:), &
                                        T1,T2,D1,D2)
  end if

  if (ANY(wantedDerivs(:,EOS_TABVT_PR) .GE. 0)) then

      call eos_tabGetSpeciesAnyTableData (species,            &
                                        speciesTemperature, &
                                        speciesDensity,     &
                                        (wantedDerivs .GE. 0),&
                                        EOS_TABULAR_P,      &
                                        maxval(wantedDerivs(:,EOS_TABVT_PR)),       &
                                        i,j,k,l,&
                                        tau, delta, &
                                        lowerBoundaryDens, &
                                        lowerBoundaryTemp,withinBoundaryTemp,withinBoundaryDens, &
                                        outData(:,EOS_TABVT_PR,:), &
                                        T1,T2,D1,D2)
  end if

  if (ANY(wantedDerivs(:,EOS_TABVT_ENTR) .GE. 0)) then

      call eos_tabGetSpeciesAnyTableData (species,            &
                                        speciesTemperature, &
                                        speciesDensity,     &
                                        (wantedDerivs .GE. 0),&
                                        EOS_TABULAR_S,      &
                                        maxval(wantedDerivs(:,EOS_TABVT_ENTR)),       &
                                        i,j,k,l,&
                                        tau, delta, &
                                        lowerBoundaryDens, &
                                        lowerBoundaryTemp,withinBoundaryTemp,withinBoundaryDens, &
                                        outData(:,EOS_TABVT_ENTR,:), &
                                        T1,T2,D1,D2)
  end if

!
!
!   ...Ready! 
!
!
  return
end subroutine eos_tabGetSpeciesTabulatedData
