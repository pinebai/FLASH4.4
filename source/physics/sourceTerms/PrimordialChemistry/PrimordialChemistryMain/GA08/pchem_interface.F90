!!****ih* source/physics/sourceTerms/PrimordialChemistry/PrimordialChemistryMain/GA08/pchem_interface
!!
!! SYNOPSIS
!!
!!  use pchem_interface
!!
!!  DESCRIPTION
!!
!!  This is the header file for the PrimordialChemistry module
!!
!!***

Module pchem_interface

#include "Flash.h"
#include "constants.h"

!!
  interface
     subroutine pchem_mapNetworkToSpecies(networkIn, specieOut)
       implicit none
       integer, intent(IN)   ::  networkIn
       integer, intent(OUT)  ::  specieOut
     end subroutine pchem_mapNetworkToSpecies
  end interface


  !!NEED BURNER HERE
  !!AND ANY OTHER SPECIFICS

  interface
     subroutine pchem_burner(tstep,temp,density,xIn,xOut,sdotRate,ei,counts,jcounts,mfrac,tback)
       implicit none
       real, intent(IN)                         :: tstep,temp,density
       real, intent(OUT)                        :: sdotRate
       real, intent(IN), dimension(NSPECIES)    :: xIn
       real, intent(OUT), dimension(NSPECIES)   :: xOut
       real, intent(INOUT)                      :: ei
       real, intent(INOUT)                      :: counts,jcounts
       real, intent(IN)                         :: mfrac
       real, intent(OUT)                        :: tback
     end subroutine pchem_burner
  end interface

  interface 
     subroutine pchem_azbar()
       implicit none
     end subroutine pchem_azbar
  end interface

  interface
     subroutine pchem_gift(ab,n1,n2)
       implicit none
       integer, INTENT(IN) :: n1,n2
       real, INTENT(INOUT) :: ab(n1,n2)
     end subroutine pchem_gift
  end interface

! going to add my new fuctions for the dense matrix solvers

  interface
     subroutine ludcmp(a,n,np,indx,d)
       implicit none
       integer                  :: n,np,indx(*)
       real                     :: a(*), d
     end subroutine ludcmp
  end interface

  interface
     subroutine lubksb(a,n,np,indx,b)
       implicit none
       integer                  :: n,np,indx(*)
       real                     :: a(*), b(*)
     end subroutine lubksb
  end interface

  interface
    subroutine  leqs(a,b,n,np)
      implicit none
      integer                   :: n,np
      real                      :: a(*),b(*)
    end subroutine leqs
  end interface


  
end Module pchem_interface
