!!****if* source/physics/Eos/EosMain/multiTemp/Helmholtz/eos_helmData
!!
!! NAME
!!
!!  eos_helmData
!!
!! 
!! SYNOPSIS
!!
!!  use eos_helmData
!!
!! DESCRIPTION
!!
!!  General parameters (non-array) for EOS Helmholtz
!!
!!
!!*** 

module eos_helmData

  implicit none

#include "constants.h"
#include "Flash.h"
#include "Eos.h"  
#include "Eos_map.h"

  ! maximum number of iterations for the Newton loop to find T from e
  integer, save :: eos_maxNewton
 
  ! general grid parameters
!!  integer, save :: eos_meshMe, eos_meshNumProcs 
  
  ! how accurately to do the Newton iteration to get T from e
  real, save :: eos_tol
  
  !real, save :: eos_smallt
  
  integer,save :: eos_hfetInit

  ! force the iterative solver to leave inputs alone (always true in MODE_DENS_TEMP)
!!MANOS  logical, save :: eos_forceConstantInput

  ! Coulomb multiplier 
  real, save :: eos_coulombMult
  ! abort if pressures become negative
  logical, save :: eos_coulombAbort
 
  integer,parameter :: EOSIMAX=211,EOSJMAX=71


  real, save :: eos_tlo, eos_tstpi
  real, save :: eos_dlo, eos_dstpi
  real,dimension(EOSJMAX),save :: eos_dt,eos_dtSqr,eos_dtInv,eos_dtSqrInv,eos_t
  real,dimension(EOSIMAX),save :: eos_dd,eos_ddSqr,eos_ddInv,eos_ddSqrInv,eos_d

!..for the helmholtz free energy tables
!..for the pressure derivative with density tables
!..for the chemical potential tables
!..for the number density tables
  real,save,dimension(EOSIMAX,EOSJMAX) :: eos_f,eos_fd, eos_ft,eos_fdd,&
                                          eos_ftt,eos_fdt,eos_fddt,&
                                          eos_fdtt, eos_fddtt, & 
                                          eos_dpdf,eos_dpdfd,eos_dpdft,&
                                          eos_dpdfdd,eos_dpdftt,eos_dpdfdt,&
                                          eos_ef,eos_efd,eos_eft,eos_efdd,&
                                          eos_eftt,eos_efdt, & 
                                          eos_xf,eos_xfd,eos_xft,eos_xfdd,&
                                          eos_xftt,eos_xfdt
end module eos_helmData
