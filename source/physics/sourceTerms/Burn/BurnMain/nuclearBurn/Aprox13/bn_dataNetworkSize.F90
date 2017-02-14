!!****if* source/physics/sourceTerms/Burn/BurnMain/nuclearBurn/Aprox13/bn_dataNetworkSize
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
  integer, parameter :: nrat = 59
  integer, parameter :: nratp1 = nrat+1

  !! For communication between bn_networkSparseJakob and bn_networkSparsePointers
  !! were in common block elca13
  integer, parameter      ::   neloc=65
  integer, save           ::   eloc(neloc),nterms

end Module bn_dataNetworkSize
