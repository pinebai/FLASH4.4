!!****if* source/physics/materialProperties/Opacity/OpacityMain/OPAL/op_opalData
!!
!! NAME
!!
!!  op_opalData
!!
!! SYNOPSIS
!!
!!  use op_opalData
!!
!! DESCRIPTION
!!
!!  Defines and stores the local data for the tabulated opacities.
!!  
!!***

#include "Flash.h"

module op_opalData
  
  implicit none

  integer, parameter :: OP_OPAL_LOWT  = 1
  integer, parameter :: OP_OPAL_HIGHT = 2

  real,    save :: op_opalMaxLowT = 1.0e4

  logical, save :: op_initializedTabulated = .false.
  logical, save :: op_useLogTables

  integer, save :: op_maxTablesLowT
  integer, save :: op_maxTablesHighT
  integer, save :: op_maxTablesRO

  character (len=80), allocatable, save :: op_tableKind      (:)
  character (len=80), allocatable, save :: op_tableNameLowT  (:)
  character (len=80), allocatable, save :: op_tableNameHighT (:)

  real,    allocatable, save :: op_opalTableAbundMax     (:)

  type opT_oneVarTablePT
     type(opT_varTableGroupPT),pointer :: pg ! group to which this table belongs
     real, pointer                      :: table(:,:) ! the data
     logical                            :: isLogData
     character(len=80)                  :: fromFile !for debugging
     integer                            :: tableNo  !for debugging
     integer                            :: lineNo   !for debugging
     integer                            :: tempRange ! whether (1) LowT, (2) HighT
     integer                            :: component ! whether (1) EOS_TAB_FOR_ION,
                                                     ! (2) EOS_TAB_FOR_ELE, (3) EOS_TAB_FOR_MAT, ?..
!     integer                            :: iSave     ! cached location of previous lookup
!     integer                            :: kSave     ! cached location of previous lookup
  end type opT_oneVarTablePT

  ! enhanced variant of opT_oneVarTablePT, for multigroup variables, think opacities.
  ! Not currently used in Eos code.
  type opT_mgVarTablePT
     type(opT_varTableGroupPT),pointer :: pg ! group to which this table belongs
     real, pointer                      :: table(:,:,:) ! the data  <- this is different from opT_oneVarTablePT
     logical                            :: isLogData
     character(len=80)                  :: fromFile !for debugging
     integer                            :: tableNo  !for debugging
     integer                            :: lineNo   !for debugging
     integer                            :: tempRange ! whether (1) LowT, (2) HighT
!     integer                            :: iSave     ! cached location of previous lookup
!     integer                            :: kSave     ! cached location of previous lookup
  end type opT_mgVarTablePT

  type opT_tableGroupDescT
     integer                           :: ntemp
     integer                           :: ndens
     integer                           :: nmg ! length of vector in multigroup variables
     real, pointer                     :: Temperatures(:)
     real, pointer                     :: Densities(:)
     logical                           :: isLog !whether *temperatures* and *densities* are stored as logarithms.
  end type opT_tableGroupDescT

  ! Now a type that can contain several 2-dimensional data tables, each of them
  ! of the same size and row and columns ranges.
  ! We shall use one object of this type for each species for each of (Z, eint, pres, hc)
  ! where necessary.
  type opT_varTableGroupPT
     type(opT_oneVarTablePT), pointer :: table ! table for one or several variable
     type(opT_mgVarTablePT), pointer  :: mgTable(:) ! tables for one or several multigroup variables
     integer                           :: numTables !DEV: maybe redundant - KW
     type(opT_tableGroupDescT)        :: td
  end type opT_varTableGroupPT

  ! Now a type that contains all the tables (or pointers to them) that
  ! pertain to one material.
  type opT_varAllTablesPT
     type(opT_varTableGroupPT)  :: tg !!! (1:OP_OPAL_NUM_H_ABUNDANCES ?)
  end type opT_varAllTablesPT

  
  type(opT_varAllTablesPT), allocatable,target,save :: op_opalAllTab(:,:) ! (OP_OPAL_NUM_H_ABUNDANCES,OP_OPAL_LOWT:OP_OPAL_HIGHT)

end module op_opalData

