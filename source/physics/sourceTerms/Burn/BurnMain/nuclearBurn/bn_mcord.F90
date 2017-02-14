!!****if* source/physics/sourceTerms/Burn/BurnMain/nuclearBurn/bn_mcord
!!
!! NAME
!!   bn_mcord
!!
!! SYNOPSIS
!!
!!   call bn_mcord(integer(IN)  ::i,
!!              integer(IN)  ::j,
!!              integer(OUT) ::iloc(:),
!!              integer(OUT) ::jloc(:),
!!              integer(INOUT)  ::nzo,
!!              integer(IN)  ::np,    DEV this is a guess on intent for nzo, np and np2, nterm
!!              integer(OUT) ::eloc(:),
!!              integer(INOUT) ::nterm,
!!              integer(IN)  ::np2)   
!!
!! DESCRIPTION
!!
!!  routine mcord counts the sparse matrix pointers
!!  marks the locations of a sparse matrix, and the contributions to 
!!  the matrix element. input is the matrix location i,j
!!  output is the is the vector iloc of physical dimension np,
!!  jloc of physical dimension np, eloc of physical dimension np2,
!!  the number of sparse matrix elements nzo, the number of contributing
!!  terms nterm.
!!  
!!  the variables nzo and nterm must be initialized (set to zero)
!!  before the first call to this routine.
!!
!! ARGUMENTS
!!
!!   i -      matrix location
!!   j -      matrix location
!!   iloc -   vector of dimension np
!!   jloc -   vector of dimension np
!!   nzo -    number of spare matrix elements
!!   np -     dimension of vector iloc, jloc
!!   eloc -   vector of dimension np2
!!   nterm -  number of contributing terms
!!   np2 -    dimension of vector eloc
!!
!! NOTES
!!   Called by the network routine bn_networkSparsePointers
!!
!!***



subroutine bn_mcord(i,j,iloc,jloc,nzo,np,eloc,nterm,np2)
  use Driver_interface, ONLY : Driver_abortFlash

  implicit none

  !!  declare arguments
  integer, intent(IN)  ::  i,j, np, np2
  integer, intent(INOUT) :: nterm, nzo
  integer, intent(OUT) ::  iloc(np),jloc(np),eloc(np2)

  !! declare local variables
  integer              :: n


  !!  increment the number of terms
  nterm = nterm + 1
  if (nterm .gt. np2) then  
     write(6,*) 
     write(6,*) 'nterm=',nterm,' np2=',np2
     write(6,*) 'nterm > np2 in routine mcord'
     write(6,*) 'check np2 dimensioning in calling routine'
     write(6,*) 
     call Driver_abortFlash('ERROR in bn_mcord: check np2 dimensioning in caller')
     stop
  end if

  !!  if the location has already been tagged, mark eloc and return
  do n=1,nzo
     if (iloc(n) .eq. i  .and.  jloc(n) .eq. j) then
        eloc(nterm) = n
        return
     end if
  enddo

  !!  location has been tagged, increment and add it to the list
  nzo          = nzo + 1
  if (nzo .gt. np) then
     write(6,*) 
     write(6,*) 'nzo=',nzo,'  np=',np
     write(6,*) 'nzo > np in routine mcord'
     write(6,*) 'check np dimensioning in calling routine'
     write(6,*) 
     call Driver_abortFlash('ERROR in bn_mcord: check np dimensioning in caller')
     stop
  end if
  eloc(nterm) = nzo
  iloc(nzo)   = i
  jloc(nzo)   = j

  return

end subroutine bn_mcord


