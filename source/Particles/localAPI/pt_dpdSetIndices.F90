!!****if* source/Particles/localAPI/pt_dpdSetIndices
!!
!! NAME
!!
!!  pt_dpdSetIndices
!!
!! SYNOPSIS
!!
!!  call pt_dpdSetIndices(integer(OUT) :: opind(3),
!!                        integer(OUT) :: npind(3),
!!                        integer(OUT) :: ovind(3),
!!                        integer(OUT) :: nvind(3),
!!                        integer(OUT) :: intvind(3),
!!                        integer(OUT) :: ofind(3),
!!                        integer(OUT) :: nfind(3))
!!
!! DESCRIPTION
!!
!! Stub 
!!
!! ARGUMENTS
!!
!!   opind : 
!!
!!   npind : 
!!
!!   ovind : 
!!
!!   nvind : 
!!
!!   intvind : 
!!
!!   ofind : 
!!
!!   nfind : 
!!
!!
!!
!!***

subroutine pt_dpdSetIndices(opind,npind,ovind,nvind,intvind,ofind,nfind)

  implicit none
#include "Flash.h"
#include "constants.h"
  
  integer,dimension(3),INTENT(out)::opind,npind,ovind,nvind,intvind,ofind,nfind
  
 ! Set the indices for variables at n and n+1
  ! Old positions and velocities indices
  opind(:)= 0
  ovind(:)= 0
  ofind(:)= 0
  ! N positions and velocities indices
  npind(:)= 0
  nvind(:)= 0
  nfind(:)= 0
  intvind(:)= 0
  

end subroutine pt_dpdSetIndices
