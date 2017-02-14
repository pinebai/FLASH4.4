!!****ih* source/physics/Eos/localAPI/eos_tabInterface
!!
!! NAME
!!     eos_tabInterface
!!
!! SYNOPSIS
!!     use eos_tabInterface
!!
!! DESCRIPTION
!!
!! This is an interface module for internal use of
!! the Tabulated implementation of the Eos unit.
!!
!!***

module eos_tabInterface
  implicit none
#include "Eos.h"
!#include "Flash.h"

  interface
     subroutine eos_tabBrowseIonmixTables (tableName,                   &
                                       needZFTable,                 &
                                       needENTable,                 &
                                       needHCTable,                 &
                                               nstepsDensityZF,     &
                                               nstepsDensityEN,     &
                                               nstepsDensityHC,     &
                                               nstepsTemperatureZF, &
                                               nstepsTemperatureEN, &
                                               nstepsTemperatureHC  )

       logical,            intent (in)  :: needZFTable
       logical,            intent (in)  :: needENTable
       logical,            intent (in)  :: needHCTable
       integer,            intent (out) :: nstepsDensityZF
       integer,            intent (out) :: nstepsDensityEN
       integer,            intent (out) :: nstepsDensityHC
       integer,            intent (out) :: nstepsTemperatureZF
       integer,            intent (out) :: nstepsTemperatureEN
       integer,            intent (out) :: nstepsTemperatureHC
       character (len=80), intent (in)  :: tableName
     end subroutine eos_tabBrowseIonmixTables
  end interface

  interface
     subroutine eos_tabBrowseTables (tableKind,                   &
                                 tableName,                   &
                                 groupName,                   &
                                 needZFTable,                 &
                                 needENTable,                 &
                                 needHCTable,                 &
                                 needEntrTable,                 &
                                         nstepsDensityZF,     &
                                         nstepsDensityEN,     &
                                         nstepsDensityHC,     &
                                         nstepsDensityEntr,     &
                                         nstepsTemperatureZF, &
                                         nstepsTemperatureEN, &
                                         nstepsTemperatureHC, &
                                         nstepsTemperatureEntr  )

       logical,            intent (in)  :: needZFTable
       logical,            intent (in)  :: needENTable
       logical,            intent (in)  :: needHCTable
       logical,            intent (in)  :: needEntrTable
       integer,            intent (out) :: nstepsDensityZF
       integer,            intent (out) :: nstepsDensityEN
       integer,            intent (out) :: nstepsDensityHC
       integer,            intent (out) :: nstepsDensityEntr
       integer,            intent (out) :: nstepsTemperatureZF
       integer,            intent (out) :: nstepsTemperatureEN
       integer,            intent (out) :: nstepsTemperatureHC
       integer,            intent (out) :: nstepsTemperatureEntr
       character (len=80), intent (in)  :: tableKind
       character (len=80), intent (in)  :: tableName
       character (len=80), intent (in)  :: groupName
     end subroutine eos_tabBrowseTables
  end interface

  interface
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

  implicit none
  integer, intent (in)  :: species
  real,    intent (in)  :: speciesTemperature
  real,    intent (in)  :: speciesDensity
  integer, intent (in)  :: selTT
  integer,intent(OUT) :: i,j,k,l
  real,intent(out) :: tau,delta
  logical,intent(OUT) :: lowerBoundaryDens
  logical,intent(OUT) :: withinBoundaryDens
  logical,intent(OUT) :: lowerBoundaryTemp
  logical,intent(OUT) :: withinBoundaryTemp
  integer,intent(INOUT) :: iPrev,kPrev
  real,intent(OUT) :: D1,D2,T1,T2
end subroutine eos_tabFindTablePos
  end interface

  interface
     subroutine eos_tabGetSpeciesAnyTableData(species,           &
                                             speciesTemperature, &
                                             speciesDensity,     &
                                             wanted, &
                                             selTT             , &
                                             needDerivs        , &
                                        i,j,k,l,&
                                        tauIn, deltaIn, &
                                        lowerBoundaryDens,&
                                        lowerBoundaryTemp,withinBoundaryTemp,withinBoundaryDens,&
                                             resultTT, &
                                        T1in,T2in,D1in,D2in)
       use eos_tabData,ONLY: EOS_TAB_NCOMP,EOS_TAB_NALLTAB
       implicit none
       integer, intent (in)  :: species
       real,    intent (in)  :: speciesTemperature
       real,    intent (in)  :: speciesDensity
       logical,            intent (in) :: wanted(EOS_TAB_NCOMP,EOS_TAB_NALLTAB)
       integer, intent (in)  :: selTT
       integer, intent (in)  :: needDerivs
  integer,intent(IN) :: i,j,k,l
  real,intent(IN) :: tauIn,deltaIn
  logical,intent(IN) :: lowerBoundaryDens
  logical,intent(IN) :: withinBoundaryDens
  logical,intent(IN) :: lowerBoundaryTemp
  logical,intent(IN) :: withinBoundaryTemp
       real,    intent (out) :: resultTT(0:,:)
  real,intent(IN) :: D1in,D2in,T1in,T2in
     end subroutine eos_tabGetSpeciesAnyTableData
  end interface

  interface
     subroutine eos_tabGetSpeciesTabulatedData (species,            &
                                                 speciesTemperature, &
                                                 speciesDensity,     &
                                                 wantedDerivs, &
                                                 outData              )
       use eos_tabData,ONLY: EOS_TAB_NCOMP,EOS_TAB_NALLTAB, &
                             EOS_TABVT_ENTR, &
                             EOS_TAB_NDERIVS
       implicit none
       integer, intent (in) :: wantedDerivs(EOS_TAB_NCOMP,EOS_TAB_NALLTAB)
       integer, intent (in) :: species
       real,    intent (in) :: speciesTemperature
       real,    intent (in) :: speciesDensity
       real,intent(OUT) :: outData(0:EOS_TAB_NDERIVS,1:EOS_TABVT_ENTR,EOS_TAB_NCOMP)
     end subroutine eos_tabGetSpeciesTabulatedData
  end interface

  interface
     subroutine eos_tabIonmix4(mode, vecLen, eosData, vecBegin, vecEnd, eosType, subtype, material, mask)
       implicit none
       integer, INTENT(in) :: mode, vecLen
       real,INTENT(inout), dimension(EOS_NUM*vecLen) :: eosData 
       integer,optional,INTENT(in) :: vecBegin,vecEnd
       integer,optional,INTENT(in) :: eosType,subtype
       integer,optional,INTENT(in) :: material
       logical, optional, INTENT(in),target,dimension(EOS_VARS+1:EOS_NUM) :: mask
     end subroutine Eos_tabIonmix4
  end interface

  interface
     subroutine eos_tabOpacplot(mode, vecLen, eosData, vecBegin, vecEnd, eosType, subtype, material, mask)
       implicit none
       integer, INTENT(in) :: mode, vecLen
       real,INTENT(inout), dimension(EOS_NUM*vecLen) :: eosData 
       integer,optional,INTENT(in) :: vecBegin,vecEnd
       integer,optional,INTENT(in) :: eosType,subtype
       integer,optional,INTENT(in) :: material
       logical, optional, INTENT(in),target,dimension(EOS_VARS+1:EOS_NUM) :: mask
     end subroutine Eos_tabOpacplot
  end interface


  interface
     subroutine eos_tabLogOutsideCounts(force)
       implicit none
       logical, INTENT(in) :: force
     end subroutine Eos_tabLogOutsideCounts
  end interface

  interface
     subroutine eos_tabReadIonmixTables (tableName,   &
                                     wanted,      &
                                     td,& 
                                     tbZF,tbEN,tbHC)
       use eos_tabData,ONLY: EOS_TAB_NCOMP,EOS_TAB_NALLTAB, &
            eosT_tableGroupDescT, &
            eosT_oneVarTablePT

       implicit none

       logical,            intent (in) :: wanted(EOS_TAB_NCOMP,EOS_TAB_NALLTAB)
       type(eosT_tableGroupDescT),intent(inout) :: td(:)
       type(eosT_oneVarTablePT),pointer,dimension(:) :: tbZF,tbEN,tbHC

       character (len=80), intent (in) :: tableName
     end subroutine eos_tabReadIonmixTables
  end interface

  interface
     subroutine eos_tabReadIonmix4Tables (tableName,   &
                                     wanted,      &
                                     td,& 
                                     tbZF,tbEN,tbPR,tbHC,tbEntr)
       use eos_tabData,ONLY: EOS_TAB_NCOMP,EOS_TAB_NALLTAB, &
            eosT_tableGroupDescT, &
            eosT_oneVarTablePT
       implicit none
       logical,            intent (in) :: wanted(EOS_TAB_NCOMP,EOS_TAB_NALLTAB)
       type(eosT_tableGroupDescT),intent(inout) :: td(:)
       type(eosT_oneVarTablePT),pointer,dimension(:) :: tbZF,tbEN,tbPR,tbHC,tbEntr

       character (len=80), intent (in) :: tableName
     end subroutine eos_tabReadIonmix4Tables
  end interface

  interface
     subroutine eos_tabReadOpacplotTables (tableName,   &
                                           groupName,   &
                                     wanted,      &
                                     td,& 
                                     tbZF,tbEN,tbPR,tbHC,tbEntr)
       use eos_tabData,ONLY: EOS_TAB_NCOMP,EOS_TAB_NALLTAB, &
            eosT_tableGroupDescT, &
            eosT_oneVarTablePT
       implicit none
       logical,            intent (in) :: wanted(EOS_TAB_NCOMP,EOS_TAB_NALLTAB)
       type(eosT_tableGroupDescT),intent(inout) :: td(:)
       type(eosT_oneVarTablePT),pointer,dimension(:) :: tbZF,tbEN,tbPR,tbHC,tbEntr

       character (len=80), intent (in) :: tableName
       character (len=80), intent (in) :: groupName
     end subroutine eos_tabReadOpacplotTables
  end interface

  interface
     subroutine eos_tabReadPropaceosTables (tableName,   &
                                     wanted,      &
                                     td,& 
                                     tbZF,tbEN,tbPR,tbHC)
       use eos_tabData,ONLY: EOS_TAB_NCOMP,EOS_TAB_NALLTAB, &
            eosT_tableGroupDescT, &
            eosT_oneVarTablePT
       implicit none
       logical,            intent (in) :: wanted(EOS_TAB_NCOMP,EOS_TAB_NALLTAB)
       type(eosT_tableGroupDescT),intent(inout) :: td(:)
       type(eosT_oneVarTablePT),pointer,dimension(:) :: tbZF,tbEN,tbPR,tbHC

       character (len=80), intent (in) :: tableName
     end subroutine eos_tabReadPropaceosTables
  end interface

  interface
     subroutine eos_tabReadTables (tableKind,   &
                               tableName,   &
                               groupName,   &
                               wanted,      &
                               td,&
                                     tbZF,tbEN,tbPR,tbHC,tbEntr)
       use eos_tabData,ONLY: EOS_TAB_NCOMP,EOS_TAB_NALLTAB, &
            eosT_tableGroupDescT, &
            eosT_oneVarTablePT

       implicit none

       logical,            intent (in) :: wanted(EOS_TAB_NCOMP,EOS_TAB_NALLTAB)
  type(eosT_tableGroupDescT),intent(inout) :: td(:)
       type(eosT_oneVarTablePT),pointer,dimension(:) :: tbZF,tbEN,tbPR,tbHC,tbEntr

       character (len=80), intent (in) :: tableKind
       character (len=80), intent (in) :: tableName
       character (len=80), intent (in) :: groupName
     end subroutine eos_tabReadTables
  end interface

  interface
     subroutine eos_tabUpdateOutsideCount(species, &
                                     upperBoundaryTemp, tempIsLog, &
                                     temp, &
                                     upperBoundaryDens, densIsLog, &
                                     dens)
       implicit none
       integer,intent(in) :: species
       logical,intent(in) :: upperBoundaryTemp, upperBoundaryDens
       logical,intent(in) :: tempIsLog, densIsLog
       real,   intent(in) :: temp, dens
     end subroutine eos_tabUpdateOutsideCount
  end interface

  interface
     subroutine eos_tabWriteSpeciesZFTable (fileUnit,species)
       integer, intent (in) :: fileUnit
       integer, intent (in) :: species
     end subroutine eos_tabWriteSpeciesZFTable
  end interface

  interface
     subroutine eos_tabWriteSpeciesENTable (fileUnit,species)
       integer, intent (in) :: fileUnit
       integer, intent (in) :: species
     end subroutine eos_tabWriteSpeciesENTable
  end interface

  interface
     subroutine eos_tabWriteSpeciesHCTable (fileUnit,species)
       integer, intent (in) :: fileUnit
       integer, intent (in) :: species
     end subroutine eos_tabWriteSpeciesHCTable
  end interface

  interface
     subroutine eos_tabWriteSpeciesPRTable (fileUnit,species)
       integer, intent (in) :: fileUnit
       integer, intent (in) :: species
     end subroutine eos_tabWriteSpeciesPRTable
  end interface

  interface
     subroutine eos_tabWriteSpeciesEntrTable (fileUnit,species)
       integer, intent (in) :: fileUnit
       integer, intent (in) :: species
     end subroutine eos_tabWriteSpeciesEntrTable
  end interface

  interface
     subroutine eos_tabWriteSpeciesTables (fileUnit,    &
                                       species,     &
                                       needZFTable, &
                                       needENTable, &
                                       needHCTable, &
                                       needPRTable, &
                                       needEntrTable)

       logical, intent (in) :: needZFTable
       logical, intent (in) :: needENTable
       logical, intent (in) :: needHCTable
       logical, intent (in) :: needPRTable
       logical, intent (in) :: needEntrTable
       integer, intent (in) :: fileUnit
       integer, intent (in) :: species
     end subroutine eos_tabWriteSpeciesTables
  end interface

  interface
     subroutine eos_tabWriteTables (fileUnit,header)
       character*(*), intent (in) :: header
       integer,       intent (in) :: fileUnit
     end subroutine eos_tabWriteTables
  end interface

#if(0)
  interface
     subroutine eos_tabWriteConstants (fileUnit,header)
       character*(*), intent (in) :: header
       integer,       intent (in) :: fileUnit
     end subroutine eos_tabWriteConstants
  end interface
#endif

  interface
     subroutine eos_tabZeroOutsideCounts()
     end subroutine Eos_tabZeroOutsideCounts
  end interface

end module eos_tabInterface
