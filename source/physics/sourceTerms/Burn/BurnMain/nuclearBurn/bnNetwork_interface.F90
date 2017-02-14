!!****ih* source/physics/sourceTerms/Burn/BurnMain/nuclearBurn/bnNetwork_interface
!!
!! NAME
!!
!!  bnNetwork_interface
!!
!!
!! SYNOPSIS
!!
!!  use bnNetwork_interface, ONLY: 
!!
!!
!! DESCRIPTION
!!
!!  Subroutine argument descriptions
!!
!!***


Module bnNetwork_interface

  interface 
     subroutine bn_initNetwork
       implicit none
     end subroutine bn_initNetwork
  end interface

  interface 
     subroutine bn_networkRates
       implicit none
     end subroutine bn_networkRates
  end interface

  interface 
     subroutine bn_networkTable
       implicit none
     end subroutine bn_networkTable
  end interface

  interface 
     subroutine bn_gift(ab,n1,n2)
       implicit none
       integer, INTENT(in) :: n1,n2
       real, INTENT(inout) :: ab(n1,n2)
     end subroutine bn_gift
  end interface

  interface 
     subroutine bn_networkScreen(y)
#include "Flash.h"
       implicit none
       real, intent(IN) ::  y(NSPECIES)
     end subroutine bn_networkScreen
  end interface

  interface
     subroutine bn_networkWeak(y)
#include "Flash.h"
       implicit none
       real, intent(IN) :: y(NSPECIES)
     end subroutine bn_networkWeak
  end interface

!----------------------------------------------------
  interface 
     subroutine bn_network(tt,y,dydt)   
       implicit none
       real, intent(IN) :: tt
       real, intent(INOUT), dimension(*)  :: y
       real, intent(OUT), dimension(*) :: dydt
     end subroutine bn_network
  end interface


! This module is the dummy argument passes. 
!  It was necessary in flash two because every network
!  called a different name e.g. aprox13.F90, aprox19.F90
  interface derivs  ! = bn_network
     subroutine derivs(tt,y,dydt)   !! == bn_network
       implicit none
       real, intent(IN) :: tt
       real, intent(INOUT), dimension(*)  :: y
       real, intent(OUT), dimension(*) :: dydt
     end subroutine derivs
  end interface


!-----------------------
!  These two routines are supposed to have the same interface.  Sadly, they
!     do not -- notice how dfdy is dimensioned.  bn_networkSparseJakob got 
!     a fake dummy argument at the end to even them out
!  No matter what I do, I can't trick a good compiler (Lahey) into thinking
!     that the two have the same interface.  Lahey figures out that dfdy
!     is one-dimensional in the case of SparseJakob and two-dimensional in the
!     case of DenseJakob.  So I give up and just hope that the twit who
!     designed these routines did his job right and it doesn't matter.
!  Back to use "external bn_networkSparseJakob, bn_networkDenseJakob, derivs"
  interface
     subroutine bn_networkSparseJakob(tt,y,dfdy,nzo,nDummy)
       implicit none
       integer, intent(IN) :: nzo, nDummy ! added to have equal numbers of arguments
       real, intent(IN)    :: tt
       real, intent(INOUT) :: y(*)
       real, intent(OUT)   :: dfdy(*)
     end subroutine bn_networkSparseJakob
  end interface

  interface 
     subroutine bn_networkDenseJakob(tt,y,dfdy,nlog,nphys)   
       implicit none
       integer, intent(IN) :: nlog, nphys
       real, intent(IN)    :: tt
       real, intent(INOUT) ::  y(*)
       real, intent(OUT)   ::  dfdy(nphys,nphys)
     end subroutine bn_networkDenseJakob
  end interface

!! Dummy version
  interface 
     subroutine jakob(tt,y,dfdy,nzo,nDummy) ! = bn_networkSparseJakob or bn_networkDenseJakob
       implicit none
       integer, intent(IN) :: nzo, nDummy
       real, intent(IN)    :: tt
       real, intent(INOUT) :: y(*)
       real, intent(OUT)   :: dfdy(nzo,nDummy)
     end subroutine jakob
  end interface

!----------------------------------------------
! This is the routine passed as bjakob. Bizarrely enough, there is no
!!  equivalent of bn_networkDensePointers
  interface
     subroutine bn_networkSparsePointers(iloc,jloc,nzo,np)
       implicit none
       integer, intent(IN)  ::   iloc(*),jloc(*),np
       integer, intent(OUT) ::   nzo
     end subroutine bn_networkSparsePointers
  end interface

! Dummy version
  interface
     subroutine bjakob(iloc,jloc,nzo,np)
       implicit none
       integer, intent(IN)  ::   iloc(*),jloc(*),np
       integer, intent(OUT) ::   nzo
     end subroutine bjakob
  end interface

end Module bnNetwork_interface

