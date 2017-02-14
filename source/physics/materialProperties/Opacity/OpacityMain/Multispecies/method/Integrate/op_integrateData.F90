!!****if* source/physics/materialProperties/Opacity/OpacityMain/Multispecies/method/Integrate/op_integrateData
!!
!! NAME
!!
!!  op_integrateData
!!
!! SYNOPSIS
!!
!!  use op_integrateData
!!
!! DESCRIPTION
!!
!!  Defines and stores the local data for the opacity integration section.
!!
!! PARAMETERS
!!
!!  op_minRombergSteps : the minimum number of Romberg iterations performed for each integration
!!  op_maxRombergSteps : the maximum number of Romberg iterations performed for each integration
!!  
!!***
module op_integrateData
  
  implicit none

  logical, save :: op_initializedIntegrate = .false.
  logical, save :: op_useRomberg
  logical, save :: op_useQuadrature
  logical, save :: op_printQuadratureData

  integer, save :: op_maxRoots
  integer, save :: op_maxWork

  real,    save :: op_ebaseRescaleExponent
  real,    save :: op_RombergAccuracy

  integer, parameter :: op_minRombergSteps = 4
  integer, parameter :: op_maxRombergSteps = 24

  real, allocatable, save :: op_AuxPolynomialA       (:)
  real, allocatable, save :: op_AuxPolynomialB       (:)
  real, allocatable, save :: op_Moments              (:)
  real, allocatable, save :: op_JmatrixDiagonals     (:)
  real, allocatable, save :: op_OrthoPolynomialA     (:)
  real, allocatable, save :: op_OrthoPolynomialB     (:)
  real, allocatable, save :: op_JmatrixOffdiagonals  (:)
  real, allocatable, save :: op_work1                (:)
  real, allocatable, save :: op_work2                (:)
  real, allocatable, save :: op_work3                (:)
  real, allocatable, save :: op_RombergIntegral      (:)
  real, allocatable, save :: op_RombergRow           (:,:)

end module op_integrateData
