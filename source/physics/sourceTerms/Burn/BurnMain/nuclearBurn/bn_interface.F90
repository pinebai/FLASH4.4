!!****ih* source/physics/sourceTerms/Burn/BurnMain/nuclearBurn/bn_interface
!!
!! SYNOPSIS
!!   use bn_interface
!!
!! DESCRIPTION
!!
!! This is the header file for the Burn module that defines its
!! private interfaces.
!!
!!***

Module bn_interface

#include "Flash.h"
#include "constants.h"

!! NOTE These next three routines are specific for each alpha-chain network
  interface
     subroutine bn_mapNetworkToSpecies(networkIn,specieOut)
       implicit none
       integer, intent(IN)  ::   networkIn
       integer, intent(OUT) :: specieOut
     end subroutine bn_mapNetworkToSpecies
  end interface

  interface
     subroutine bn_burner(tstep,temp,density,xIn,xOut,sdotRate)
       implicit none
       real, intent(IN)               :: tstep,temp,density
       real, intent(OUT)              :: sdotRate
       real, intent(IN), dimension(NSPECIES) :: xIn   
       real, intent(OUT), dimension(NSPECIES):: xOut
     end subroutine bn_burner
  end interface

  interface
     subroutine bn_gift(ab,n1,n2)
       implicit none
       integer, INTENT(in) :: n1,n2
       real, INTENT(inout) :: ab(n1,n2)
     end subroutine bn_gift
  end interface

  !! NOTE from this point down, LBR had to guess the intent
  !! These routines reside in nuclearBurn and are generic for all alpha-chain networks
  interface
     subroutine bn_azbar()
       implicit none
     end subroutine bn_azbar
  end interface

  interface
     subroutine bn_ecapnuc(etakep,temp,rpen,rnep,spen,snep)
       implicit none
       real, intent(IN)     :: etakep, temp
       real, intent(OUT)    :: rpen,rnep,spen,snep
     end subroutine bn_ecapnuc
  end interface

  interface
     real function bn_ifermi12(f)
       implicit none
       real, intent(IN) :: f
     end function bn_ifermi12
  end interface

  interface
     subroutine bn_mazurek(btemp,bden,y56,ye,rn56ec,sn56ec)
       implicit none
       real, intent(IN) :: btemp, bden
       real, intent(IN) :: y56, ye
       real, intent(OUT):: rn56ec,sn56ec

     end subroutine bn_mazurek
  end interface

  interface
     subroutine bn_mcord(i,j,iloc,jloc,nzo,np,eloc,nterm,np2)
       implicit none
       integer, intent(IN)  ::  i,j, np, np2
       integer, intent(INOUT) :: nterm, nzo
       integer, intent(OUT) ::  iloc(np),jloc(np),eloc(np2)
     end subroutine bn_mcord
  end interface

  interface
     subroutine bn_screen4(zbarr,abarr,z2barr,z1,a1,z2,a2, & 
          &                   jscreen,init,scorr,scorrdt)
       implicit none
       integer, intent(IN)   :: jscreen, init
       real, intent(IN)      :: abarr, zbarr, z2barr, z1, a1, z2, a2
       real, intent(OUT)     :: scorr
       real, intent(OUT), optional     :: scorrdt
     end subroutine bn_screen4
  end interface

  interface
     subroutine bn_sneutx()
       implicit none
     end subroutine bn_sneutx
  end interface

end Module bn_interface
