!!****ih* source/physics/sourceTerms/Ionize/IonizeMain/ion_interface
!!
!! NAME
!!
!!  ion_interface
!!
!! SYNOPSIS
!!
!!  use ion_interface
!!
!! DESCRIPTION
!!
!!  Interface module for internal use within Ionize.
!!
!!***

Module ion_interface 
  implicit none

#include "constants.h"
#include "Flash.h"
#include "Ionize.h"
  
  interface
     subroutine ion_wint(xaux,n,x,j)
       implicit none
       integer, intent(IN) :: n
       real, dimension(n), intent(IN) :: xaux
       integer, intent(OUT) :: j
       real, intent(IN) :: x
     end subroutine ion_wint
  end interface

  interface
     subroutine ion_selectElements()
       implicit none
     end subroutine ion_selectElements
  end interface

  interface
     subroutine ion_readTable()
       implicit none
     end subroutine ion_readTable
  end interface

  interface
     subroutine ion_intCoeff(nel,tx,al,cz,nion)
       implicit none
       integer, intent(IN) :: nel
       real, intent(INOUT) :: tx
       real, dimension(ION_NIMAX), intent(OUT) :: al
       real, dimension(0:ION_NIMAX), intent(OUT) :: cz
       integer, intent(OUT) :: nion
     end subroutine ion_intCoeff
  end interface

  interface
     subroutine ion_mcord(i,j,iloc,jloc,nzo,np,eloc,nterm,np2)
       implicit none
       integer, intent(INOUT) :: nzo, nterm
       integer, intent(IN)    :: i, j, np, np2
       integer, dimension(np),  intent(INOUT)  :: iloc, jloc
       integer, dimension(np2), intent(INOUT)  :: eloc
     end subroutine ion_mcord
  end interface

end Module ion_interface
