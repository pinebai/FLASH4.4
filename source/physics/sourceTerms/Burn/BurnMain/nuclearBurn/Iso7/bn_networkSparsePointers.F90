!!****ih* source/physics/sourceTerms/Burn/BurnMain/nuclearBurn/Iso7/bn_networkSparsePointers
!!
!! NAME
!!  
!!  bn_networkSparsePointers
!!
!! SYNOPSIS
!! 
!!  call bn_networkSparsePointers(integer, intent(INOUT)  :: iloc(np),
!!                                integer, intent(INOUT)  :: jloc(np),
!!                                integer, intent(OUT)    :: nzo,
!!                                integer, intent(IN)     :: np)
!!
!!  
!! DESCRIPTION
!!
!!  routine networkSparsePointers builds the nonzero locations for networkSparseJakob
!!
!!  input is the integer arrys iloc and jloc, both of dimension np, that
!!  on output contain nzo matrix element locations.
!!  
!!  This routine is called as an external from bn_burner, as routine 'bjakob'
!!   
!!  ARGUMENTS
!!  
!!  iloc   -- integer array of size np
!!  jloc   -- integer array of size np  
!!  np     -- integer, size of arrays
!!  nzo    -- integer, number of nonzero matrix element locations
!!
!!***

subroutine bn_networkSparsePointers(iloc,jloc,nzo,np)

  use Burn_data
  use bn_dataIso7
  use bn_dataNetworkSize, ONLY : neloc, eloc, nterms 

  implicit none
  !      save

#include "Flash.h"

  !.. 
  !..this routine builds the nonzero matrix locations for bn_bn_networkSparseJakob
  !..input is the integer arrys iloc and jloc, both of dimension np, that
  !..on output contain nzo matrix element locations.
  !.. 

  !..declare arguments
  integer, intent(IN)  ::   iloc(*),jloc(*),np
  integer, intent(OUT) ::   nzo

  !..declare locals
  integer          :: i


  !..communicate with bn_bn_networkSparseJakob in bn_dataIso7


  !..initialize
  nterms = 0
  nzo    = 0
  do i=1,neloc 
     eloc(i) = 0
  enddo

  !..tag the nonzero locations
  !..4he jacobian elementss
  call bn_mcord(ihe4,ihe4,iloc,jloc,nzo,np,eloc,nterms,neloc)
  call bn_mcord(ihe4,ic12,iloc,jloc,nzo,np,eloc,nterms,neloc)
  call bn_mcord(ihe4,io16,iloc,jloc,nzo,np,eloc,nterms,neloc)
  call bn_mcord(ihe4,ine20,iloc,jloc,nzo,np,eloc,nterms,neloc)
  call bn_mcord(ihe4,img24,iloc,jloc,nzo,np,eloc,nterms,neloc)
  call bn_mcord(ihe4,isi28,iloc,jloc,nzo,np,eloc,nterms,neloc)
  call bn_mcord(ihe4,ini56,iloc,jloc,nzo,np,eloc,nterms,neloc)

  !..12c jacobian elements 
  call bn_mcord(ic12,ihe4,iloc,jloc,nzo,np,eloc,nterms,neloc)
  call bn_mcord(ic12,ic12,iloc,jloc,nzo,np,eloc,nterms,neloc)
  call bn_mcord(ic12,io16,iloc,jloc,nzo,np,eloc,nterms,neloc)

  !..16o jacobian elements 
  call bn_mcord(io16,ihe4,iloc,jloc,nzo,np,eloc,nterms,neloc)
  call bn_mcord(io16,ic12,iloc,jloc,nzo,np,eloc,nterms,neloc)
  call bn_mcord(io16,io16,iloc,jloc,nzo,np,eloc,nterms,neloc)
  call bn_mcord(io16,ine20,iloc,jloc,nzo,np,eloc,nterms,neloc)

  !..20ne jacobian elements 
  call bn_mcord(ine20,ihe4,iloc,jloc,nzo,np,eloc,nterms,neloc)
  call bn_mcord(ine20,ic12,iloc,jloc,nzo,np,eloc,nterms,neloc)
  call bn_mcord(ine20,io16,iloc,jloc,nzo,np,eloc,nterms,neloc)
  call bn_mcord(ine20,ine20,iloc,jloc,nzo,np,eloc,nterms,neloc)
  call bn_mcord(ine20,img24,iloc,jloc,nzo,np,eloc,nterms,neloc)

  !..24mg jacobian elements 
  call bn_mcord(img24,ihe4,iloc,jloc,nzo,np,eloc,nterms,neloc)
  call bn_mcord(img24,ic12,iloc,jloc,nzo,np,eloc,nterms,neloc)
  call bn_mcord(img24,io16,iloc,jloc,nzo,np,eloc,nterms,neloc)
  call bn_mcord(img24,ine20,iloc,jloc,nzo,np,eloc,nterms,neloc)
  call bn_mcord(img24,img24,iloc,jloc,nzo,np,eloc,nterms,neloc)
  call bn_mcord(img24,isi28,iloc,jloc,nzo,np,eloc,nterms,neloc)

  !..28si jacobian elements 
  call bn_mcord(isi28,ihe4,iloc,jloc,nzo,np,eloc,nterms,neloc)
  call bn_mcord(isi28,ic12,iloc,jloc,nzo,np,eloc,nterms,neloc)
  call bn_mcord(isi28,io16,iloc,jloc,nzo,np,eloc,nterms,neloc)
  call bn_mcord(isi28,img24,iloc,jloc,nzo,np,eloc,nterms,neloc)
  call bn_mcord(isi28,isi28,iloc,jloc,nzo,np,eloc,nterms,neloc)
  call bn_mcord(isi28,ini56,iloc,jloc,nzo,np,eloc,nterms,neloc)

  !..ni56 jacobian elements 
  call bn_mcord(ini56,ihe4,iloc,jloc,nzo,np,eloc,nterms,neloc)
  call bn_mcord(ini56,isi28,iloc,jloc,nzo,np,eloc,nterms,neloc)
  call bn_mcord(ini56,ini56,iloc,jloc,nzo,np,eloc,nterms,neloc)


  !..write a diagnostic
  !..      write(6,*) ' '
  !..      write(6,*) nzo,' matrix elements'
  !..      write(6,*) nterms,' jacobian contributions'
  !..      write(6,*) ' '

  return
end subroutine bn_networkSparsePointers


