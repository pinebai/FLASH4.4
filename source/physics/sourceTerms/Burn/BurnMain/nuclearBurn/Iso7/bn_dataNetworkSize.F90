!!****if* source/physics/sourceTerms/Burn/BurnMain/nuclearBurn/Iso7/bn_dataNetworkSize
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
!!    nrat = number of reaction rates
!!
!! NOTES
!!
!!   In Flash2, this routine was called network_size.fh
!!
!!***

Module bn_dataNetworkSize

  integer          nrat,nratp1
  parameter        (nrat = 15, nratp1 = nrat+1)

!! For communication between bn_networkSparePointers and bn_networkSparseJakob
!! were in common block elcpp
!!  common /elcpp/  eloc,nterms
  integer, parameter      ::   neloc=34
  integer, save           ::   eloc(neloc),nterms


end Module bn_dataNetworkSize
