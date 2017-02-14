!!****ih* source/physics/sourceTerms/Burn/BurnMain/nuclearBurn/Aprox13/bn_networkSparsePointers
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
  ! for communication with networkSparseJakob
  use bn_dataNetworkSize, ONLY : neloc, eloc, nterms

  use bn_interface, ONLY: bn_mcord

  implicit none

  !    save

  !!  declare
  integer, intent(INOUT) ::   iloc(*),jloc(*)
  integer, intent(IN)    ::   np
  integer, intent(OUT)   ::   nzo

  !!  local variables
  integer                ::    i


  !!  communicate with networkSparseJakob, now in bn_dataAprox13
!!  integer          neloc
!!  parameter        (neloc=65)
!!  integer          eloc(neloc),nterms
!!  common /elca13/  eloc,nterms


  !!  initialize
  nterms = 0
  nzo    = 0
  do i=1,neloc 
     eloc(i) = 0
  enddo

  !!  tag the nonzero locations
  !!  he4 jacobian elements
  call bn_mcord(ihe4,ihe4,iloc,jloc,nzo,np,eloc,nterms,neloc)
  call bn_mcord(ihe4,ic12,iloc,jloc,nzo,np,eloc,nterms,neloc)
  call bn_mcord(ihe4,io16,iloc,jloc,nzo,np,eloc,nterms,neloc)
  call bn_mcord(ihe4,ine20,iloc,jloc,nzo,np,eloc,nterms,neloc)
  call bn_mcord(ihe4,img24,iloc,jloc,nzo,np,eloc,nterms,neloc)
  call bn_mcord(ihe4,isi28,iloc,jloc,nzo,np,eloc,nterms,neloc)
  call bn_mcord(ihe4,is32,iloc,jloc,nzo,np,eloc,nterms,neloc)
  call bn_mcord(ihe4,iar36,iloc,jloc,nzo,np,eloc,nterms,neloc)
  call bn_mcord(ihe4,ica40,iloc,jloc,nzo,np,eloc,nterms,neloc)
  call bn_mcord(ihe4,iti44,iloc,jloc,nzo,np,eloc,nterms,neloc)
  call bn_mcord(ihe4,icr48,iloc,jloc,nzo,np,eloc,nterms,neloc)
  call bn_mcord(ihe4,ife52,iloc,jloc,nzo,np,eloc,nterms,neloc)
  call bn_mcord(ihe4,ini56,iloc,jloc,nzo,np,eloc,nterms,neloc)

  !!  c12 jacobian elements
  call bn_mcord(ic12,ihe4,iloc,jloc,nzo,np,eloc,nterms,neloc)
  call bn_mcord(ic12,ic12,iloc,jloc,nzo,np,eloc,nterms,neloc)
  call bn_mcord(ic12,io16,iloc,jloc,nzo,np,eloc,nterms,neloc)

  !!  16o jacobian elements
  call bn_mcord(io16,ihe4,iloc,jloc,nzo,np,eloc,nterms,neloc)
  call bn_mcord(io16,ic12,iloc,jloc,nzo,np,eloc,nterms,neloc)
  call bn_mcord(io16,io16,iloc,jloc,nzo,np,eloc,nterms,neloc)
  call bn_mcord(io16,ine20,iloc,jloc,nzo,np,eloc,nterms,neloc)

  !!  20ne jacobian elements
  call bn_mcord(ine20,ihe4,iloc,jloc,nzo,np,eloc,nterms,neloc)
  call bn_mcord(ine20,ic12,iloc,jloc,nzo,np,eloc,nterms,neloc)
  call bn_mcord(ine20,io16,iloc,jloc,nzo,np,eloc,nterms,neloc)
  call bn_mcord(ine20,ine20,iloc,jloc,nzo,np,eloc,nterms,neloc)
  call bn_mcord(ine20,img24,iloc,jloc,nzo,np,eloc,nterms,neloc)

  !!  24mg jacobian elements
  call bn_mcord(img24,ihe4,iloc,jloc,nzo,np,eloc,nterms,neloc)
  call bn_mcord(img24,ic12,iloc,jloc,nzo,np,eloc,nterms,neloc)
  call bn_mcord(img24,io16,iloc,jloc,nzo,np,eloc,nterms,neloc)
  call bn_mcord(img24,ine20,iloc,jloc,nzo,np,eloc,nterms,neloc)
  call bn_mcord(img24,img24,iloc,jloc,nzo,np,eloc,nterms,neloc)
  call bn_mcord(img24,isi28,iloc,jloc,nzo,np,eloc,nterms,neloc)

  !!  28si jacobian elements
  call bn_mcord(isi28,ihe4,iloc,jloc,nzo,np,eloc,nterms,neloc)
  call bn_mcord(isi28,ic12,iloc,jloc,nzo,np,eloc,nterms,neloc)
  call bn_mcord(isi28,io16,iloc,jloc,nzo,np,eloc,nterms,neloc)
  call bn_mcord(isi28,img24,iloc,jloc,nzo,np,eloc,nterms,neloc)
  call bn_mcord(isi28,isi28,iloc,jloc,nzo,np,eloc,nterms,neloc)
  call bn_mcord(isi28,is32,iloc,jloc,nzo,np,eloc,nterms,neloc)

  !!  32s jacobian elements
  call bn_mcord(is32,ihe4,iloc,jloc,nzo,np,eloc,nterms,neloc)
  call bn_mcord(is32,io16,iloc,jloc,nzo,np,eloc,nterms,neloc)
  call bn_mcord(is32,isi28,iloc,jloc,nzo,np,eloc,nterms,neloc)
  call bn_mcord(is32,is32,iloc,jloc,nzo,np,eloc,nterms,neloc)
  call bn_mcord(is32,iar36,iloc,jloc,nzo,np,eloc,nterms,neloc)

  !!  36ar jacobian elements
  call bn_mcord(iar36,ihe4,iloc,jloc,nzo,np,eloc,nterms,neloc)
  call bn_mcord(iar36,is32,iloc,jloc,nzo,np,eloc,nterms,neloc)
  call bn_mcord(iar36,iar36,iloc,jloc,nzo,np,eloc,nterms,neloc)
  call bn_mcord(iar36,ica40,iloc,jloc,nzo,np,eloc,nterms,neloc)

  !!  40ca jacobian elements
  call bn_mcord(ica40,ihe4,iloc,jloc,nzo,np,eloc,nterms,neloc)
  call bn_mcord(ica40,iar36,iloc,jloc,nzo,np,eloc,nterms,neloc)
  call bn_mcord(ica40,ica40,iloc,jloc,nzo,np,eloc,nterms,neloc)
  call bn_mcord(ica40,iti44,iloc,jloc,nzo,np,eloc,nterms,neloc)

  !!  44ti jacobian elements
  call bn_mcord(iti44,ihe4,iloc,jloc,nzo,np,eloc,nterms,neloc)
  call bn_mcord(iti44,ica40,iloc,jloc,nzo,np,eloc,nterms,neloc)
  call bn_mcord(iti44,iti44,iloc,jloc,nzo,np,eloc,nterms,neloc)
  call bn_mcord(iti44,icr48,iloc,jloc,nzo,np,eloc,nterms,neloc)

  !!  48cr jacobian elements
  call bn_mcord(icr48,ihe4,iloc,jloc,nzo,np,eloc,nterms,neloc)
  call bn_mcord(icr48,iti44,iloc,jloc,nzo,np,eloc,nterms,neloc)
  call bn_mcord(icr48,icr48,iloc,jloc,nzo,np,eloc,nterms,neloc)
  call bn_mcord(icr48,ife52,iloc,jloc,nzo,np,eloc,nterms,neloc)

  !!  52fe jacobian elements
  call bn_mcord(ife52,ihe4,iloc,jloc,nzo,np,eloc,nterms,neloc)
  call bn_mcord(ife52,icr48,iloc,jloc,nzo,np,eloc,nterms,neloc)
  call bn_mcord(ife52,ife52,iloc,jloc,nzo,np,eloc,nterms,neloc)
  call bn_mcord(ife52,ini56,iloc,jloc,nzo,np,eloc,nterms,neloc)

  !!  56ni jacobian elements
  call bn_mcord(ini56,ihe4,iloc,jloc,nzo,np,eloc,nterms,neloc)
  call bn_mcord(ini56,ife52,iloc,jloc,nzo,np,eloc,nterms,neloc)
  call bn_mcord(ini56,ini56,iloc,jloc,nzo,np,eloc,nterms,neloc)


  !!  
  !!  write a diagnostic
  !      write(6,*) ' '
  !      write(6,*) nzo,' matrix elements'
  !      write(6,*) nterms,' jacobian contributions'
  !      write(6,*) ' '

  return

end   subroutine bn_networkSparsePointers
