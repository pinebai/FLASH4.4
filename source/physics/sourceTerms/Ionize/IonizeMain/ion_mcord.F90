!!****if* source/physics/sourceTerms/Ionize/IonizeMain/ion_mcord
!!
!!  NAME
!!
!!   ion_mcord
!!
!! SYNOPSIS
!!
!!   ion_mcord(integer (IN)    :: i,
!!             integer (IN)    :: j,
!!             integer (INOUT) :: iloc(np),
!!             integer (INOUT) :: jloc(np),
!!             integer (INOUT) :: nzo,
!!             integer (IN)    :: np,
!!             integer (INOUT) :: eloc(np2),
!!             integer (INOUT) :: nterm,
!!             integer (IN)    :: np2)
!!
!!  DESCRIPTION
!!
!!    marks the locations of a sparse matrix, and the contributions to
!!    the matrix element. input is the matrix location i,j
!!    output is the vector iloc of physical dimension np,
!!    jloc of physical dimension np, eloc of physical dimension np2,
!!    the number of sparse matrix elements nzo, the number of contributing
!!    terms nterm.
!! 
!! ARGUMENTS
!!  
!!   i : i index of matrix location
!!   j : j index of matrix location
!!   iloc(:) : vector of tagged indices i
!!   jloc(:) : vector of tagged indices j
!!   nzo : number of matrix elements
!!   np :  physical dimension of iloc and jloc
!!   eloc(:) : number of matrix elements for a given nterm
!!   nterm : number of contributing terms
!!   np2 : physical dimension for eloc
!!
!!
!! NOTES
!!   the variables nzo and nterm must be initialized (set to zero)
!!   before the first call to this routine.
!!***

subroutine ion_mcord(i,j,iloc,jloc,nzo,np,eloc,nterm,np2)
    
  use Driver_interface, ONLY: Driver_abortFlash

  implicit none
    

  !..declare
  integer, intent(INOUT) :: nzo, nterm
  integer, intent(IN)    :: i, j, np, np2
  integer, dimension(np),  intent(INOUT)  :: iloc, jloc
  integer, dimension(np2), intent(INOUT)  :: eloc
  integer                :: n
  !
  !..increment the number of terms
  nterm = nterm + 1
  if (nterm .gt. np2) then
     write(6,*)
     write(6,*) 'nterm=',nterm,' np2=',np2
     write(6,*) 'nterm > np2 in routine mcord'
     write(6,*) 'check np2 dimensioning in calling routine'
     write(6,*)
     call Driver_abortFlash('Error in mcord: check np2 dimmension in calling routine')
  end if
  !
  !..if the location has already been tagged, mark eloc and return
  do n = 1,nzo
     if (iloc(n) .eq. i  .and.  jloc(n) .eq. j) then
        eloc(nterm) = n
        return
     end if
  enddo
  !
  !..location has been tagged, increment and add it to the list
  nzo = nzo + 1
  if (nzo .gt. np) then
     write(6,*)
     write(6,*) 'nzo=',nzo,'  np=',np
     write(6,*) 'nzo > np in routine mcord'
     write(6,*) 'check np dimensioning in calling routine'
     write(6,*)
     call Driver_abortFlash('Error in mcord: check mp dimensioning in calling routine')
  end if
  eloc(nterm) = nzo
  iloc(nzo)   = i
  jloc(nzo)   = j
  return
end subroutine ion_mcord
