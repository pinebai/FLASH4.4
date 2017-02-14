!!****ih* source/physics/sourceTerms/PrimordialChemistry/PrimordialChemistryMain/GA08/pchem_networkInterface
!!
!! NAME
!!
!!  chemNetwork_interface
!!
!! SYNOPSIS
!!
!!   use chemNetwork_interface
!!
!! DESCRIPTION
!!
!!  Subroutine argument descriptions
!!
!!***


Module pchem_networkInterface

  interface
     subroutine pchem_initNetwork
       implicit none
     end subroutine pchem_initNetwork
  end interface


  interface
     subroutine pchem_networkRates
       implicit none
     end subroutine pchem_networkRates
  end interface


  interface
     subroutine pchem_networkScreen(y)
#include "Flash.h"
       implicit none
       real, intent(IN) :: y(NSPECIES)
     end subroutine pchem_networkScreen
  end interface
!!ABOVE MAY NOT WORK. y(*) should be y(NSPECIES)


  interface
     subroutine pchem_network(tt,y,dydt)
       implicit none
       real, intent(IN) :: tt
       real, intent(INOUT), dimension(*) :: y
       real, intent(OUT), dimension(*) :: dydt
     end subroutine pchem_network
  end interface

  interface
     subroutine pchem_gift(ab,n1,n2)
       implicit none
       integer, INTENT(IN)  :: n1,n2
       real, INTENT(INOUT)  :: ab(n1,n2)
     end subroutine pchem_gift
  end interface


! Don't know if I need this but will put it in. 
! From burn

  interface derivs
     subroutine derivs(tt,y,dydt)
       implicit none
       real, intent(IN) :: tt
     !  real, intent(INOUT), dimension(*) :: y
     !  real, intent(OUT), dimension(*)   :: dydt
       real, intent(INOUT) :: y(*)
       real, intent(OUT)   :: dydt(*)
     end subroutine derivs
  end interface


  interface
     subroutine pchem_networkDenseJakob(tt,y,dfdy,nlog,nphys)
       implicit none
       integer, intent(IN) :: nlog, nphys
       real, intent(IN) :: tt
       real, intent(INOUT) :: y(*)
     !  real, intent(OUT)   :: dfdy(nlog,nlog)
       real, intent(OUT)   :: dfdy(nphys,nphys)
     end subroutine pchem_networkDenseJakob
  end interface


  interface
     subroutine jakob(tt,y,dfdy,nzo,nDummy)
       implicit none
       integer, intent(IN)  :: nzo, nDummy
       real, intent(IN)     :: tt
       real, intent(INOUT)  :: y(*)
       real, intent(OUT)    :: dfdy(nzo,nDummy)
    end subroutine jakob
  end interface
  
  
  interface 
     subroutine pchem_networkSparsePointers(iloc,jloc,nzo,np)
       implicit none
       integer, intent(IN)  :: iloc(*),jloc(*),np
       integer, intent(OUT) :: nzo
     end subroutine pchem_networkSparsePointers
  end interface


  interface
     subroutine bjakob(iloc,jloc,nzo,np)
       implicit none
       integer, intent(IN)  :: iloc(*),jloc(*),np
       integer, intent(OUT) :: nzo
     end subroutine
  end interface

end Module pchem_networkInterface