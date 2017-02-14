!!****if* source/physics/materialProperties/Opacity/OpacityMain/Multispecies/method/LowTemp/op_lowTempData
!!
!! NAME
!!
!!  op_lowTempData
!!
!! SYNOPSIS
!!
!!  use op_lowTempData
!!
!! DESCRIPTION
!!
!!  Defines and stores the local data for the low temperature opacity section.
!!
!! PARAMETERS
!!
!!  op_maxElements      : the maximum number of atomic elements that can be handled
!!                        (from atomic numbers Z = 1 to op_maxElements)
!!  op_maxNstepsLowTemp : the size of the low temperature grid
!!  
!!***
module op_lowTempData
  
  implicit none

  logical, save :: op_ignoreKleinNishina
  logical, save :: op_ignorePhotoElectron
  logical, save :: op_initializedLowTemp = .false.

  integer, save :: op_maxJmax
  real,    save :: op_maxLowTempTemperature
  real,    save :: op_minLowTempTemperature

  integer, parameter :: op_maxElements      = 100
  integer, parameter :: op_maxNstepsLowTemp =  30

  real,     allocatable, save :: op_A1group                (:)
  real,     allocatable, save :: op_A2group                (:)
  real,     allocatable, save :: op_A3group                (:)
  real,     allocatable, save :: op_A4group                (:)
  real,     allocatable, save :: op_intLimits              (:)

  integer,  allocatable, save :: op_Jmax                   (:)
  real,     allocatable, save :: op_PEenergyRange          (:,:,:)
  real,     allocatable, save :: op_Aij4                   (:,:,:)

  integer,  allocatable, save :: op_elementJmax            (:)
  real,     allocatable, save :: op_elementPEenergyRange   (:,:,:)
  real,     allocatable, save :: op_elementAij4            (:,:,:)

  real,     allocatable, save :: op_tableLowTemp           (:)
  real,     allocatable, save :: op_PlanckLowTempTables    (:,:,:)
  real,     allocatable, save :: op_RosselandLowTempTables (:,:,:)

end module op_lowTempData
