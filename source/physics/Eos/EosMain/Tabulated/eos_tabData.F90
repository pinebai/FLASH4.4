!!****if* source/physics/Eos/EosMain/Tabulated/eos_tabData
!!
!! NAME
!!
!!  eos_tabData
!!
!! SYNOPSIS
!!  use eos_tabData
!!
!! DESCRIPTION
!!
!!  Defines and stores the local data for Unit Opacity.
!!  
!!***

#include "Eos.h"
#include "Eos_components.h"

module eos_tabData
  
  implicit none

#if (EOSCOMP_NUM_COMPONENTS==1)
  integer, parameter :: EOS_TAB_FOR_ION = 1
  integer, parameter :: EOS_TAB_FOR_ELE = 1
  integer, parameter :: EOS_TAB_FOR_MAT = 1
#else
  integer, parameter :: EOS_TAB_FOR_ION = 1
  integer, parameter :: EOS_TAB_FOR_ELE = 2
  integer, parameter :: EOS_TAB_FOR_MAT = 3
#endif

  integer, parameter :: EOS_TAB_NCOMP = EOS_TAB_FOR_MAT

  ! variable types - families of tabulated variables
  integer, parameter :: EOS_TABVT_ZF = 1 !Z_free
  integer, parameter :: EOS_TABVT_EN = 2 !(internal) energies
  integer, parameter :: EOS_TABVT_PR = 3 !pressures
  integer, parameter :: EOS_TABVT_ENTR = 4 !entropies
  integer, parameter :: EOS_TABVT_HC = 5 !heat capacities
  integer, parameter :: EOS_TABVT_OP = 6 !opacities
  integer, parameter :: EOS_TAB_NALLTAB = 6

  integer, parameter :: EOS_TABINT_DERIV_0 = 0 ! leave as 0
  integer, parameter :: EOS_TABINT_DERIV_DT = 1 ! leave consecutive
  integer, parameter :: EOS_TABINT_DERIV_DD = 2 ! leave consecutive
  integer, parameter :: EOS_TAB_NDERIVS = 2


  logical, save :: eos_useLogTables

  integer, save :: op_maxNstepsDensityZF
  integer, save :: op_maxNstepsDensityEN
  integer, save :: op_maxNstepsDensityHC
  integer, save :: op_maxNstepsTemperatureZF
  integer, save :: op_maxNstepsTemperatureEN
  integer, save :: op_maxNstepsTemperatureHC
  integer, save :: op_maxTablesZF
  integer, save :: op_maxTablesEN
  integer, save :: op_maxTablesHC
  integer, save :: eos_tabTotalNumSpecies


  real, parameter :: zero =  0.0
  real, parameter :: one  =  1.0
  real, parameter :: ten  = 10.0

  character (len=80), allocatable, save ::  eos_tableKind (:)
  character (len=80), allocatable, save ::  eos_tableName (:)
  character (len=80), allocatable, save ::  eos_groupName (:)

  integer, allocatable, save :: eos_tabIonizationKind    (:,:)
  integer, allocatable, save :: eos_tabIntEnergyKind     (:,:)
  integer, allocatable, save :: eos_tabHeatCpKind        (:,:)


  type eosT_oneVarTablePT
     type(eosT_varTableGroupPT),pointer :: pg ! group to which this table belongs
     real, pointer                      :: table(:,:) ! the data
     logical                            :: isLogData
     character(len=80)                  :: fromFile !for debugging
     integer                            :: tableNo  !for debugging
     integer                            :: lineNo   !for debugging
     integer                            :: derivedFrom1 !for debugging ?
     integer                            :: derivedFrom2 !for debugging ?
     integer                            :: varType ! whether (1) Z, (2) eint, (3) pres, (4) hc, ..
     integer                            :: component ! whether (1) EOS_TAB_FOR_ION,
                                                     ! (2) EOS_TAB_FOR_ELE, (3) EOS_TAB_FOR_MAT, ?..
!     integer                            :: iSave     ! cached location of previous lookup
!     integer                            :: kSave     ! cached location of previous lookup
  end type eosT_oneVarTablePT

  ! enhanced variant of eosT_oneVarTablePT, for multigroup variables, think opacities.
  ! Not currently used in Eos code.
  type eosT_mgVarTablePT
     type(eosT_varTableGroupPT),pointer :: pg ! group to which this table belongs
     real, pointer                      :: table(:,:,:) ! the data  <- this is different from eosT_oneVarTablePT
     logical                            :: isLogData
     character(len=80)                  :: fromFile !for debugging
     integer                            :: tableNo  !for debugging
     integer                            :: lineNo   !for debugging
     integer                            :: derivedFrom1 !for debugging ?
     integer                            :: derivedFrom2 !for debugging ?
     integer                            :: varType ! whether absoption / emission / Rosseland ?
     integer                            :: component ! whether (1) for ions,
                                                     ! (2) for electrons, (3) for matter?
!     integer                            :: iSave     ! cached location of previous lookup
!     integer                            :: kSave     ! cached location of previous lookup
  end type eosT_mgVarTablePT

  type eosT_tableGroupDescT
     integer                           :: ntemp
     integer                           :: ndens
     integer                           :: nmg ! length of vector in multigroup variables
     real, pointer                     :: Temperatures(:)
     real, pointer                     :: Densities(:)
     logical                           :: isLog !whether *temperatures* and *densities* are stored as logarithms.
  end type eosT_tableGroupDescT

  ! Now a type that can contain several 2-dimensional data tables, each of them
  ! of the same size an row and columns ranges.
  ! We shall use one object of this type for each species for each of (Z, eint, pres, hc)
  ! where necessary.
  type eosT_varTableGroupPT
     type(eosT_oneVarTablePT), pointer :: table(:) ! tables for one or several variables
     type(eosT_mgVarTablePT), pointer  :: mgTable(:) ! tables for one or several multigroup variables
     integer                           :: numTables !DEV: maybe redundant - KW
     type(eosT_tableGroupDescT)        :: td
!     integer                           :: ntemp
!     integer                           :: ndens
!     integer                           :: nmg ! length of vector in multigroup variables
!     real, pointer                     :: Temperatures(:)
!     real, pointer                     :: Densities(:)
  end type eosT_varTableGroupPT

  ! Now a type that contains all the tables (or pointers to them) that
  ! pertain to one material.
  type eosT_varAllTablesPT
     type(eosT_varTableGroupPT)  :: tg(1:EOS_TAB_NALLTAB)
  end type eosT_varAllTablesPT

  
  type(eosT_varAllTablesPT), allocatable,target,save :: eos_allTab(:)

  type eosT_diagEvent
     real :: temp
     real :: dens
  end type eosT_diagEvent

  type eosT_diagPT
     integer :: highTempCount
     integer :: highDensCount
     real :: highestTemp
     real :: highestDens
     type(eosT_diagEvent) :: firstHighTempEvent
     type(eosT_diagEvent) :: firstHighDensEvent
     logical :: highTempVarsLookedUp(1:EOS_TAB_NALLTAB)
     logical :: highDensVarsLookedUp(1:EOS_TAB_NALLTAB)
  end type eosT_diagPT


  type(eosT_diagPT), allocatable,target,save :: eos_tabAllDiag(:)



end module eos_tabData

