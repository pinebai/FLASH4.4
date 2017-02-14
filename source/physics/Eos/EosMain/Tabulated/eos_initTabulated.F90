!!****if* source/physics/Eos/EosMain/Tabulated/eos_initTabulated
!!
!! NAME
!!
!!  eos_initTabulated
!!
!! SYNOPSIS
!!
!!  call eos_initTabulated()
!!                      
!!                    
!!
!! DESCRIPTION
!!
!!  Initializes data for Tabulated implementations of the Eos unit.
!!
!! ARGUMENTS
!!
!!  
!!  
!!  
!!
!!***

#include "Eos.h"
#include "Flash.h"
#include "Multispecies.h"
#include "constants.h"

subroutine eos_initTabulated()

  use eos_tabData,                ONLY : EOS_TAB_FOR_ION,           &
                                         EOS_TAB_FOR_ELE,           &
                                         EOS_TAB_FOR_MAT,           &
                                         EOS_TAB_NDERIVS,           &
                                         EOS_TAB_NALLTAB,           &
                                         EOS_TABVT_ZF   ,           &
                                         EOS_TABVT_EN   ,           &
                                         EOS_TABVT_HC   ,           &
                                         EOS_TABVT_PR   ,           &
                                         EOS_TABVT_ENTR ,           &
                                          op_maxNstepsDensityZF,     &
                                          op_maxNstepsDensityEN,     &
                                          op_maxNstepsDensityHC,     &
                                          op_maxNstepsTemperatureZF, &
                                          op_maxNstepsTemperatureEN, &
                                          op_maxNstepsTemperatureHC, &
                                          op_maxTablesZF,            &
                                          op_maxTablesEN,            &
                                          op_maxTablesHC,            &
                                          EOS_TAB_NCOMP,          &
                                          eos_tabTotalNumSpecies,           &
                                          eos_useLogTables,           &
                                          eos_tableKind,              &
                                          eos_tableName,              &
                                          eos_groupName,              &
                                          eos_tabIonizationKind,         &
                                          eos_tabIntEnergyKind,           &
                                          eos_tabHeatCpKind,          &
                                          eos_allTab, &
                                          eos_tabAllDiag

  use Multispecies_interface,      ONLY : Multispecies_getProperty, Multispecies_setProperty, &
                                          Multispecies_list
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
  use PhysicalConstants_interface, ONLY : PhysicalConstants_get
  use Driver_interface,            ONLY : Driver_abortFlash
  use eos_tabInterface,            ONLY : eos_tabBrowseTables, eos_tabReadTables, &
                                          eos_tabWriteTables
  
  use Eos_data,                    ONLY : eos_globalMe, eos_meshMe, eos_entrEleScaleChoice

  implicit none

  character*80 :: dummyLine
  character*80 :: printoutHeader

  logical :: needZFTable
  logical :: needENTable
  logical :: needHCTable
  logical :: needPRTable
  logical :: needEntrTable
  logical :: needTable
  logical :: isIonmix4Like
  logical :: isIonmix4, isIonmix6
  logical :: isOpacplot
  logical :: wanted(EOS_TAB_NCOMP,EOS_TAB_NALLTAB)

  integer :: i

  integer :: fileUnit
  integer :: nstepsDensityZF
  integer :: nstepsDensityEN
  integer :: nstepsDensityHC
  integer :: nstepsDensityEntr
  integer :: nstepsTemperatureZF
  integer :: nstepsTemperatureEN
  integer :: nstepsTemperatureHC
  integer :: nstepsTemperatureEntr
  integer :: ut_getFreeFileUnit
  integer :: specno, spec

  integer :: eosType(NSPECIES), eosSubType(NSPECIES)

  real,pointer :: temperatures(:)
  real,pointer :: densities(:)
!
!
!    ...Set internal parameters.
!
!

  eos_tabTotalNumSpecies  = max (1,NSPECIES)
!
!
!    ...Allocate those arrays that depend on the # of species and # of energy groups alone.
!
!
  allocate (eos_tableKind             (1:eos_tabTotalNumSpecies))
  allocate (eos_tableName             (1:eos_tabTotalNumSpecies))
  allocate (eos_groupName             (1:eos_tabTotalNumSpecies))
  allocate (eos_tabIonizationKind    (EOS_TAB_NCOMP,1:eos_tabTotalNumSpecies))
  allocate (eos_tabIntEnergyKind     (EOS_TAB_NCOMP,1:eos_tabTotalNumSpecies))
  allocate (eos_tabHeatCpKind        (EOS_TAB_NCOMP,1:eos_tabTotalNumSpecies))

!
!
!    ...Get the external data.
!
!
  call RuntimeParameters_get ("eos_useLogTables",   eos_useLogTables )
  call RuntimeParameters_get("eos_entrEleScaleChoice", eos_entrEleScaleChoice)


  ! LOOP # 0 - Allocation & zeroing of pointers

  allocate(eos_allTab(NSPECIES))
  do specno = 1,NSPECIES
     do i = 1,EOS_TAB_NALLTAB
        !nullify pointers for tables of all types
        nullify(eos_allTab(specno)%tg(i)%table)
        nullify(eos_allTab(specno)%tg(i)%mgTable)
        eos_allTab(specno)%tg(i)%numTables = 0 ! to change seen below
        eos_allTab(specno)%tg(i)%td%ntemp = 0 ! to change seen below
        eos_allTab(specno)%tg(i)%td%ndens = 0 ! to change seen below
        eos_allTab(specno)%tg(i)%td%nmg = 0
        nullify(eos_allTab(specno)%tg(i)%td%Temperatures)
        nullify(eos_allTab(specno)%tg(i)%td%Densities)
        if (i <= EOS_TABVT_PR .OR. i==EOS_TABVT_ENTR) then
           eos_allTab(specno)%tg(i)%td%isLog = eos_useLogTables
        else
           eos_allTab(specno)%tg(i)%td%isLog = .FALSE.
        end if
     end do
  end do
!
!
!    ...Initialize all data to zero.
!
!
  op_maxTablesZF            = 0
  op_maxTablesEN            = 0
  op_maxTablesHC            = 0
  op_maxNstepsDensityZF     = 0
  op_maxNstepsDensityEN     = 0
  op_maxNstepsDensityHC     = 0
  op_maxNstepsTemperatureZF = 0
  op_maxNstepsTemperatureEN = 0
  op_maxNstepsTemperatureHC = 0

  eos_tabIonizationKind         = 0
  eos_tabIntEnergyKind           = 0
  eos_tabHeatCpKind          = 0

  wanted = .FALSE.
!
!
!    ...Assume that Multispecies_init as initialized the relevant fields in
!    the multispecies database.
!    !!DEV: Some of that initialization should be moved to the Eos unit -
!    filenames, at least. - KW
!
!


!
!
!    ...Read in the needed info from the external file and set the handle
!       arrays.
!
!
  ! LOOP # 1 - Get info from the Multispecies database that
  ! used to be read from the EOS_sources file.

  do specno = 1,eos_tabTotalNumSpecies

     spec = SPECIES_BEGIN - 1 + specno

     call Multispecies_getProperty (spec , MS_EOSTYPE    , eosType(specno))
     call Multispecies_getProperty (spec , MS_EOSSUBTYPE , eosSubType(specno))
     if (eosType(specno)==EOS_TAB) then
        select case (eosSubType(specno))
        case(4)
           eos_tableKind(specno) = "IONMIX4"
        case(6)
           eos_tableKind(specno) = "IONMIX6"
        case(1)
           eos_tableKind(specno) = "IONMIX"
        case(7)
           eos_tableKind(specno) = "OPACPLOT"
        case default
           eos_tableKind(specno) = "???"
        end select
     else if (eosType(specno)==EOS_GAM) then
           eos_tableKind(specno) = "???"
     else
           eos_tableKind(specno) = "???"
     end if
     
     isOpacplot  = (eosType(specno)==EOS_TAB .AND. eosSubType(specno)==7)
     isIonmix4   = (eosType(specno)==EOS_TAB .AND. eosSubType(specno)==4)
     isIonmix6   = (eosType(specno)==EOS_TAB .AND. eosSubType(specno)==6)

     if (eosType(specno)==EOS_TAB) then
        eos_tabIonizationKind (:,specno) = EOS_TABULAR_Z
        eos_tabIntEnergyKind (:,specno) = EOS_TABULAR_E
        if (eosSubType(specno)==4 .or.eosSubType(specno)==6) then
           eos_tabHeatCpKind (:,specno) = EOS_TABULAR_P
        else
           eos_tabHeatCpKind (:,specno) = EOS_TABULAR_C
        end if
     else
        eos_tabIonizationKind (:,specno) = EOS_APPROX_KIN
        eos_tabIntEnergyKind (:,specno) = EOS_APPROX_KIN
        eos_tabHeatCpKind (:,specno) = EOS_APPROX_KIN
     end if

     call Multispecies_getProperty (spec , MS_EOSPRESFILE , eos_tableName(specno))
     call Multispecies_getProperty (spec , MS_EOSGROUPNAME , eos_groupName(specno))

  end do
!
!
  ! LOOP # 2 - Browse & Allocate
  !
!    ...Determine first the overall maximal dimensions needed for allocating
!       the tables.
!
!    ...Allocate the tables and associated data for interpolation.
!
  do specno = 1,eos_tabTotalNumSpecies

     isOpacplot  = (eos_tableKind(specno)=='OPACPLOT')
     isIonmix4   = (eos_tableKind(specno)=='IONMIX4')
     isIonmix6   = (eos_tableKind(specno)=='IONMIX6')
     isIonmix4Like = (isIonmix4 .OR. isIonmix6)

     needZFTable  =  ANY((eos_tabIonizationKind (:,specno) == EOS_TABULAR_Z) &
                    .or. (eos_tabIntEnergyKind   (:,specno) == EOS_TABULAR_Z) &
                    .or. (eos_tabHeatCpKind  (:,specno) == EOS_TABULAR_Z) )

     needENTable  =  ANY((eos_tabIonizationKind (:,specno) == EOS_TABULAR_E) &
                    .or. (eos_tabIntEnergyKind   (:,specno) == EOS_TABULAR_E) &
                    .or. (eos_tabHeatCpKind  (:,specno) == EOS_TABULAR_E))

     needHCTable  =  ANY((eos_tabIonizationKind (:,specno) == EOS_TABULAR_C) &
                    .or. (eos_tabIntEnergyKind   (:,specno) == EOS_TABULAR_C) &
                    .or. (eos_tabHeatCpKind  (:,specno) == EOS_TABULAR_C))!!!   .AND. .FALSE. !!!!!!

     needPRTable  =  ANY((eos_tabIonizationKind (:,specno) == EOS_TABULAR_P) &
                    .or. (eos_tabIntEnergyKind   (:,specno) == EOS_TABULAR_P) &
                    .or. (eos_tabHeatCpKind  (:,specno) == EOS_TABULAR_P))

     needEntrTable = isIonmix6

     needTable    =       needZFTable &
                    .or.  needENTable &
                    .or.  needHCTable &
                    .or.  needPRTable &
                    .or.  needEntrTable

     if (needZFTable) op_maxTablesZF = op_maxTablesZF + 1
     if (needENTable) op_maxTablesEN = op_maxTablesEN + 1
     if (needHCTable) op_maxTablesHC = op_maxTablesHC + 1
     if (needPRTable) op_maxTablesHC = op_maxTablesHC + 1

     if (needTable) then

         call eos_tabBrowseTables (eos_tableKind (specno),      &
                               eos_tableName (specno),      &
                               eos_groupName (specno),      &
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
                                       nstepsTemperatureEntr)

         op_maxNstepsDensityZF     = max (nstepsDensityZF     , op_maxNstepsDensityZF    )
         op_maxNstepsDensityEN     = max (nstepsDensityEN     , op_maxNstepsDensityEN    )
         op_maxNstepsDensityHC     = max (nstepsDensityHC     , op_maxNstepsDensityHC    )
         op_maxNstepsTemperatureZF = max (nstepsTemperatureZF , op_maxNstepsTemperatureZF)
         op_maxNstepsTemperatureEN = max (nstepsTemperatureEN , op_maxNstepsTemperatureEN)
         op_maxNstepsTemperatureHC = max (nstepsTemperatureHC , op_maxNstepsTemperatureHC)

         if (needZFTable) then
            eos_allTab(specno)%tg(EOS_TABVT_ZF)%td%ntemp = nstepsTemperatureZF
            eos_allTab(specno)%tg(EOS_TABVT_ZF)%td%ndens = nstepsDensityZF
            nullify(eos_allTab(specno)%tg(EOS_TABVT_ZF)%td%Temperatures)
            nullify(eos_allTab(specno)%tg(EOS_TABVT_ZF)%td%Densities)
!!$            allocate(eos_allTab(specno)%tg(EOS_TABVT_ZF)%td%Temperatures(nstepsTemperatureZF))
!!$            allocate(eos_allTab(specno)%tg(EOS_TABVT_ZF)%td%Densities(nstepsDensityZF))
            allocate(eos_allTab(specno)%tg(EOS_TABVT_ZF)%table(EOS_TAB_FOR_ELE:EOS_TAB_FOR_ELE))
            eos_allTab(specno)%tg(EOS_TABVT_ZF)%numTables = 1
            do i = LBOUND(eos_allTab(specno)%tg(EOS_TABVT_ZF)%table,1), &
                 UBOUND(eos_allTab(specno)%tg(EOS_TABVT_ZF)%table,1)
               nullify(eos_allTab(specno)%tg(EOS_TABVT_ZF)%table(i)%table)
               eos_allTab(specno)%tg(EOS_TABVT_ZF)%table(i)%isLogData = .FALSE.
            end do
         End if
         if (needENTable) then
            eos_allTab(specno)%tg(EOS_TABVT_EN)%td%ntemp = nstepsTemperatureEN
            eos_allTab(specno)%tg(EOS_TABVT_EN)%td%ndens = nstepsDensityEN
            nullify(eos_allTab(specno)%tg(EOS_TABVT_EN)%td%Temperatures)
            nullify(eos_allTab(specno)%tg(EOS_TABVT_EN)%td%Densities)
!!$            allocate(eos_allTab(specno)%tg(EOS_TABVT_EN)%td%Temperatures(nstepsTemperatureEN))
!!$            allocate(eos_allTab(specno)%tg(EOS_TABVT_EN)%td%Densities(nstepsDensityEN))
            if (isIonmix4Like) then
               allocate(eos_allTab(specno)%tg(EOS_TABVT_EN)%table(EOS_TAB_NCOMP))
               eos_allTab(specno)%tg(EOS_TABVT_EN)%numTables = EOS_TAB_NCOMP
            else
               allocate(eos_allTab(specno)%tg(EOS_TABVT_EN)%table(EOS_TAB_FOR_MAT:EOS_TAB_FOR_MAT))
               eos_allTab(specno)%tg(EOS_TABVT_EN)%numTables = 1 
            end if
            do i = LBOUND(eos_allTab(specno)%tg(EOS_TABVT_EN)%table,1), &
                 UBOUND(eos_allTab(specno)%tg(EOS_TABVT_EN)%table,1)
               nullify(eos_allTab(specno)%tg(EOS_TABVT_EN)%table(i)%table)
               eos_allTab(specno)%tg(EOS_TABVT_EN)%table(i)%isLogData = .FALSE.
            end do
         end if
         if (needHCTable) then
            eos_allTab(specno)%tg(EOS_TABVT_HC)%td%ntemp = nstepsTemperatureHC
            eos_allTab(specno)%tg(EOS_TABVT_HC)%td%ndens = nstepsDensityHC
            nullify(eos_allTab(specno)%tg(EOS_TABVT_HC)%td%Temperatures)
            nullify(eos_allTab(specno)%tg(EOS_TABVT_HC)%td%Densities)
!!$            allocate(eos_allTab(specno)%tg(EOS_TABVT_HC)%td%Temperatures(nstepsTemperatureHC))
!!$            allocate(eos_allTab(specno)%tg(EOS_TABVT_HC)%td%Densities(nstepsDensityHC))
            allocate(eos_allTab(specno)%tg(EOS_TABVT_HC)%table(EOS_TAB_NCOMP))
            eos_allTab(specno)%tg(EOS_TABVT_HC)%numTables = EOS_TAB_NCOMP
            do i = LBOUND(eos_allTab(specno)%tg(EOS_TABVT_HC)%table,1), &
                 UBOUND(eos_allTab(specno)%tg(EOS_TABVT_HC)%table,1)
               nullify(eos_allTab(specno)%tg(EOS_TABVT_HC)%table(i)%table)
               eos_allTab(specno)%tg(EOS_TABVT_HC)%table(i)%isLogData = .FALSE.
            end do
         end if
         if (needPRTable) then
            eos_allTab(specno)%tg(EOS_TABVT_PR)%td%ntemp = nstepsTemperatureHC !DEV: ??
            eos_allTab(specno)%tg(EOS_TABVT_PR)%td%ndens = nstepsDensityHC !DEV: ??
            nullify(eos_allTab(specno)%tg(EOS_TABVT_PR)%td%Temperatures)
            nullify(eos_allTab(specno)%tg(EOS_TABVT_PR)%td%Densities)
!!$            allocate(eos_allTab(specno)%tg(EOS_TABVT_PR)%td%Temperatures(nstepsTemperatureHC))
!!$            allocate(eos_allTab(specno)%tg(EOS_TABVT_PR)%td%Densities(nstepsDensityHC))

!!$            allocate(eos_allTab(specno)%tg(EOS_TABVT_PR)%table(EOS_TAB_FOR_ION:EOS_TAB_FOR_ELE))
            if (isIonmix4Like) then
               allocate(eos_allTab(specno)%tg(EOS_TABVT_PR)%table(EOS_TAB_NCOMP))
               eos_allTab(specno)%tg(EOS_TABVT_PR)%numTables = EOS_TAB_NCOMP
            else
               allocate(eos_allTab(specno)%tg(EOS_TABVT_PR)%table(EOS_TAB_FOR_MAT:EOS_TAB_FOR_MAT))
               eos_allTab(specno)%tg(EOS_TABVT_PR)%numTables = 1 
            end if
!!$            eos_allTab(specno)%tg(EOS_TABVT_PR)%numTables = EOS_TAB_FOR_ELE - EOS_TAB_FOR_ION + 1
            do i = LBOUND(eos_allTab(specno)%tg(EOS_TABVT_PR)%table,1), &
                 UBOUND(eos_allTab(specno)%tg(EOS_TABVT_PR)%table,1)
               nullify(eos_allTab(specno)%tg(EOS_TABVT_PR)%table(i)%table)
               eos_allTab(specno)%tg(EOS_TABVT_PR)%table(i)%isLogData = .FALSE.
            end do
         end if
         if (needEntrTable) then
            eos_allTab(specno)%tg(EOS_TABVT_ENTR)%td%ntemp = nstepsTemperatureEntr
            eos_allTab(specno)%tg(EOS_TABVT_ENTR)%td%ndens = nstepsDensityEntr
            nullify(eos_allTab(specno)%tg(EOS_TABVT_ENTR)%td%Temperatures)
            nullify(eos_allTab(specno)%tg(EOS_TABVT_ENTR)%td%Densities)
!!$            allocate(eos_allTab(specno)%tg(EOS_TABVT_ENTR)%td%Temperatures(nstepsTemperatureEntr))
!!$            allocate(eos_allTab(specno)%tg(EOS_TABVT_ENTR)%td%Densities(nstepsDensityEntr))

!!$            allocate(eos_allTab(specno)%tg(EOS_TABVT_ENTR)%table(EOS_TAB_FOR_ELE))
            allocate(eos_allTab(specno)%tg(EOS_TABVT_ENTR)%table(EOS_TAB_FOR_ELE:EOS_TAB_FOR_ELE))
            eos_allTab(specno)%tg(EOS_TABVT_ENTR)%numTables = 1
            do i = LBOUND(eos_allTab(specno)%tg(EOS_TABVT_ENTR)%table,1), &
                 UBOUND(eos_allTab(specno)%tg(EOS_TABVT_ENTR)%table,1)
               nullify(eos_allTab(specno)%tg(EOS_TABVT_ENTR)%table(i)%table)
               eos_allTab(specno)%tg(EOS_TABVT_ENTR)%table(i)%isLogData = .FALSE.
            end do
         end if

     end if

  end do
!
!
!    ...Close the 'EOS_sources.txt' file.
!
!
!
!
  ! PRE-LOOP # 3  (some checks)
  !
!
!
  if (op_maxTablesZF > 0) then

      if (op_maxNstepsDensityZF == 0) then
          call Driver_abortFlash ('[eos_initTabulated] ERROR: no ZF table density grid')
      end if
      if (op_maxNstepsTemperatureZF == 0) then
          call Driver_abortFlash ('[eos_initTabulated] ERROR: no ZF table temperature grid')
      end if

  end if

  if (op_maxTablesEN > 0) then

      if (op_maxNstepsDensityEN == 0) then
          call Driver_abortFlash ('[eos_initTabulated] ERROR: no EN table density grid')
      end if
      if (op_maxNstepsTemperatureEN == 0) then
          call Driver_abortFlash ('[eos_initTabulated] ERROR: no EN table temperature grid')
      end if

  end if

  if (op_maxTablesHC > 0) then

      if (op_maxNstepsDensityHC == 0) then
          call Driver_abortFlash ('[eos_initTabulated] ERROR: no HC table density grid')
      end if
      if (op_maxNstepsTemperatureHC == 0) then
          call Driver_abortFlash ('[eos_initTabulated] ERROR: no HC table temperature grid')
      end if

  end if

!
!
!    ...Read the tables and store the needed opacity values and all associated data
!       into the arrays.
!
!

  ! LOOP # 3 - Actually call eos_tabReadTables with appropriate arguments
  ! to get the table data that will be needed into memory.
  !
  do specno = 1,eos_tabTotalNumSpecies

     wanted = .FALSE.
     isOpacplot = (eos_tableKind(specno)=='OPACPLOT')
     isIonmix4 = (eos_tableKind(specno)=='IONMIX4')
     isIonmix6 = (eos_tableKind(specno)=='IONMIX6')
     isIonmix4Like = (isIonmix4 .OR. isIonmix6)

     needZFTable  =  ANY((eos_tabIonizationKind (:,specno) == EOS_TABULAR_Z) &
                    .or. (eos_tabIntEnergyKind   (:,specno) == EOS_TABULAR_Z) &
                    .or. (eos_tabHeatCpKind  (:,specno) == EOS_TABULAR_Z))

     needENTable  =  ANY((eos_tabIonizationKind (:,specno) == EOS_TABULAR_E) &
                    .or. (eos_tabIntEnergyKind   (:,specno) == EOS_TABULAR_E) &
                    .or. (eos_tabHeatCpKind  (:,specno) == EOS_TABULAR_E))

     needHCTable  =  ANY((eos_tabIonizationKind (:,specno) == EOS_TABULAR_C) &
                    .or. (eos_tabIntEnergyKind   (:,specno) == EOS_TABULAR_C) &
                    .or. (eos_tabHeatCpKind  (:,specno) == EOS_TABULAR_C))!!!   .AND. .FALSE. !!!!!!

     needPRTable  =  ANY((eos_tabIonizationKind (:,specno) == EOS_TABULAR_P) &
                    .or. (eos_tabIntEnergyKind   (:,specno) == EOS_TABULAR_P) &
                    .or. (eos_tabHeatCpKind  (:,specno) == EOS_TABULAR_P))

     needEntrTable = isIonmix6

     needTable    =       needZFTable &
                    .or.  needENTable &
                    .or.  needHCTable &
                    .or.  needPRTable &
                    .or.  needEntrTable

     if (needZFTable) then
         wanted(EOS_TAB_FOR_ELE,EOS_TABVT_ZF) = .TRUE.
     end if

     if (needENTable) then
         if (isIonmix4Like) then
            wanted(EOS_TAB_FOR_ION:EOS_TAB_FOR_ELE,EOS_TABVT_EN) = .TRUE.
         else
            wanted(EOS_TAB_FOR_MAT,EOS_TABVT_EN) = .TRUE.
         end if
     end if

     if (needHCTable .OR. needPRTable) then
     end if
     if (needHCTable) then
        if (isIonmix4Like) then
           wanted(EOS_TAB_FOR_ION:EOS_TAB_FOR_ELE,EOS_TABVT_HC) = .TRUE.
        else
           wanted(EOS_TAB_FOR_MAT,EOS_TABVT_HC) = .TRUE.
        end if
     end if
     if (needPRTable) then
        if (isIonmix4Like) then
           wanted(EOS_TAB_FOR_ION:EOS_TAB_FOR_ELE,EOS_TABVT_PR) = .TRUE.
        else
           wanted(EOS_TAB_FOR_MAT,EOS_TABVT_PR) = .TRUE.
        end if
!!$         wanted(EOS_TAB_FOR_ELE,EOS_TABVT_PR) = .TRUE.
     end if

     if (needEntrTable) then
        wanted(EOS_TAB_FOR_ELE,EOS_TABVT_ENTR) = .TRUE.
     end if


#if(0)
     ! Currently done in Multispecies_init - KW
     if (needZFTable) &
          call Multispecies_setProperty (SPECIES_BEGIN - 1 + specno , MS_EOSZFREEFILE , eos_tableName(specno))
     if (needENTable) &
          call Multispecies_setProperty (SPECIES_BEGIN - 1 + specno , MS_EOSENERFILE , eos_tableName(specno))
     if (needPRTable) &
          call Multispecies_setProperty (SPECIES_BEGIN - 1 + specno , MS_EOSPRESFILE , eos_tableName(specno))
     if (isIonmix4) then
        call Multispecies_setProperty (SPECIES_BEGIN - 1 + specno , MS_EOSSUBTYPE , 4)
     else if (isOpacplot) then
        call Multispecies_setProperty (SPECIES_BEGIN - 1 + specno , MS_EOSSUBTYPE , 7)
     if (isIonmix6) then
        call Multispecies_setProperty (SPECIES_BEGIN - 1 + specno , MS_EOSSUBTYPE , 6)
     else
        call Multispecies_setProperty (SPECIES_BEGIN - 1 + specno , MS_EOSSUBTYPE , 1)
     end if
#endif


     nullify(temperatures)
     nullify(densities)

     if (needTable) then

        if(eos_meshMe < 4) then
           print *, "in eos_inittabulated, tableName = ", trim(eos_tableName(specno))
           print *, "in eos_inittabulated, groupName = ", trim(eos_groupName(specno))
        end if

        call eos_tabReadTables (eos_tableKind (specno), &
                                eos_tableName (specno), &
                                eos_groupName (specno), &
                                wanted,                 &
                                eos_allTab(specno)%tg(:)%td, &
                                eos_allTab(specno)%tg(EOS_TABVT_ZF)%table, &
                                eos_allTab(specno)%tg(EOS_TABVT_EN)%table, &
                                eos_allTab(specno)%tg(EOS_TABVT_PR)%table, & 
                                eos_allTab(specno)%tg(EOS_TABVT_HC)%table, &
                                eos_allTab(specno)%tg(EOS_TABVT_ENTR)%table)
     end if

  end do
!
!
!    DEV: ?? ...Print out the EOS data constants to see what has been stored.
!
!
#if(0)
  fileUnit = ut_getFreeFileUnit ()

  open (unit = fileUnit, &
        file = "opacity_printout_constants.txt", &
        form = 'formatted')

  printoutHeader = "   EOS_APPROX_KIN CONST PRINTOUT (in ?some? units)"

  call eos_tabWriteConstants (fileUnit,printoutHeader)

  close (fileUnit)
#endif
!
!
!    ...Print out the data tables to see what has been stored.
!
!

  if ( eos_globalMe == MASTER_PE ) then

     fileUnit = ut_getFreeFileUnit ()

     open (unit = fileUnit, &
          file = "EOS_printout_tables.txt", &
          form = 'formatted')

     if (eos_useLogTables) then
        printoutHeader = "   EOS TABLES PRINTOUT (logarithmic base 10, energies in units of erg/g)"
     else
        printoutHeader = "   EOS TABLES PRINTOUT (energies in units of erg/g)"
     end if

     call eos_tabWriteTables (fileUnit,printoutHeader)

     write(fileUnit,FMT='(/4x,"Multispecies_list REPRISE"/)')
     call Multispecies_list(fileUnit)

     close (fileUnit)

  end if

  ! LOOP # 4 - Initialize for EOS table diagnostics.
  !
  allocate(eos_tabAllDiag(NSPECIES))
  do specno = 1,NSPECIES
     eos_tabAllDiag(specno)%highTempCount = 0
     eos_tabAllDiag(specno)%highDensCount = 0
     eos_tabAllDiag(specno)%highestTemp = -999.0
     eos_tabAllDiag(specno)%highestDens = -999.0
     eos_tabAllDiag(specno)%highTempVarsLookedUp(:) = .FALSE.
     eos_tabAllDiag(specno)%highDensVarsLookedUp(:) = .FALSE.
     eos_tabAllDiag(specno)%firstHighTempEvent%temp = -999.0
     eos_tabAllDiag(specno)%firstHighTempEvent%dens = -999.0
     eos_tabAllDiag(specno)%firstHighDensEvent%temp = -999.0
     eos_tabAllDiag(specno)%firstHighDensEvent%dens = -999.0
  end do
  
!
!
!    ...Ready!
!
!
  return
end subroutine eos_initTabulated
