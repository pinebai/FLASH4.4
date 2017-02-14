!!****if* source/physics/sourceTerms/Burn/BurnMain/nuclearBurn/Aprox19/bn_dataNetworkSize
!!
!! NAME
!!  
!!  bn_dataNetworkSize
!!
!!
!! SYNOPSIS
!! 
!!  use bn_dataNetworkSize
!!
!! DESCRIPTION
!!
!!  Contains variables indicating the nuclear network size.
!!
!! NOTES
!!
!!   In Flash2, this routine was called network_size.fh
!!
!!***

Module bn_dataNetworkSize
  implicit none
  integer, parameter :: nrat = 86
  integer, parameter :: nratp1 = nrat+1

!! For communication between bn_networkSparseJakob and bn_networkSparsePonters
!!  were in the common /elccm1/  eloc,nterms
  integer, parameter ::     neloc=112
  integer          eloc(neloc),nterms



end Module bn_dataNetworkSize
