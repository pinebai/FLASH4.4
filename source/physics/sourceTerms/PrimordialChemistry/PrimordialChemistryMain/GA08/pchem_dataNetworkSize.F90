!!****if* source/physics/sourceTerms/PrimordialChemistry/PrimordialChemistryMain/GA08/pchem_dataNetworkSize
!!
!! NAME
!!
!! pchem_dataNetworkSize
!!
!! DESCRIPTION
!!
!! Contains variables indicating the chemistry network size
!!
!! NOTES
!!
!!  Contains variables that pertain to the size of the PrimordialChemistry Network
!!  Specifically the number of reactions
!!
!!***

Module pchem_dataNetworkSize
 
  implicit none
  integer, parameter :: nrat = 125 !Number of reaction rates
  integer, parameter :: nratp1 = nrat + 1

!! There are a couple of other terms in Burn that I think I will
!! need. But I am not sure what they are defining. 
  integer, parameter :: neloc = 65
  integer, save     :: eloc(neloc), nterms

end Module pchem_dataNetworkSize