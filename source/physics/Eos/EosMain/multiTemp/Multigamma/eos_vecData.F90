!!****if* source/physics/Eos/EosMain/multiTemp/Multigamma/eos_vecData
!!
!! NAME
!!
!!  eos_vecData
!! 
!! SYNOPSIS
!!
!!  use eos_vecData
!!
!! DESCRIPTION
!!
!!   Contains arrays for EOS Helmholtz computation.
!!   Incorporates both fixed block size and non-fixed blocksize computations
!!
!! ARGUMENTS
!!
!!
!!*** 

module eos_vecData

#include "Flash.h"


#ifdef FIXEDBLOCKSIZE
  integer, parameter :: NROWMAX = GRID_IHI_GC*GRID_JHI_GC*GRID_KHI_GC
  !..thermodynamic and composition inputs
  real, dimension(NROWMAX), save ::  tempRow,denRow,                    &
       &                             abarRow,zbarRow
  real, dimension(NROWMAX), save ::  tempRadRow,tempIonRow,tempEleRow

  !..totals and their derivatives
  real, dimension(NROWMAX), save ::  ptotRow,dptRow,             &
       &                             etotRow,detRow, stotRow            
  real, dimension(NROWMAX), save ::  dedRow, dstRow,dsdRow,dpdRow            
  real, dimension(NROWMAX), save ::  deaRow, dezRow  !Calhoun            

  !..electron-positron contributions -- most UNUSED and REMOVED
  real, dimension(NROWMAX), save :: pelRow, neRow, etaRow

  !..derivative based quantities
  real, dimension(NROWMAX), save ::    gamcRow
  real, dimension(NROWMAX), save ::    cpRow,cvRow 

#include "Eos.h"
#include "Eos_components.h"

  real, dimension(EOSCOMP_NUM_COMPONENTS,NROWMAX), save ::    eCompRow, pcompRow

#else
  !! NONFIXEDBLOCKSIZE begins here
  !
  !..thermodynamic and composition inputs
  real, allocatable, dimension(:), save ::  tempRow,denRow,                    &
       &                             abarRow,zbarRow
  real, allocatable, dimension(:), save ::  tempRadRow,tempIonRow,tempEleRow

  !..totals and their derivatives
  real, allocatable, dimension(:), save ::  ptotRow,dptRow,             &
       &                             etotRow,detRow, stotRow
  real, allocatable, dimension(:), save ::  dsdRow, dstRow, dedRow,dpdRow
  real, allocatable, dimension(:), save ::  deaRow, dezRow ! DL following Calhoun

  !..electron-positron contributions -- most UNUSED and REMOVED from allocation
  real, allocatable, dimension(:), save :: pelRow, neRow, etaRow

  !..derivative based quantities
  real, allocatable, dimension(:), save ::   gamcRow
  real, allocatable, dimension(:), save ::    cpRow,cvRow

  real, allocatable, dimension(:,:), save ::    eCompRow, pcompRow

#endif

end module eos_vecData
