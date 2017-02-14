!!****if* source/Particles/ParticlesMain/unitTest/Particles_unitTest
!!
!! NAME
!!
!!  Particles_unitTest
!!
!! SYNOPSIS
!!
!!  call Particles_unitTest(
!!                          integer(in)     :: fileUnit,
!!                          logical(inout)  :: perfect  )
!!
!! DESCRIPTION
!!
!!  This implementation tests whether particles have moved outside the local bounding
!!  box.  The simulation unitTest/Particles generates random particles and
!!  moves them about.  When this routine is called, they should all be in
!!  "sensible" places.
!!
!! ARGUMENTS
!!
!!  
!!  fileUnit:  integer   number of file for output
!!  perfect:   logical   remains true if all tests are passed
!!
!! NOTES
!!
!!***

subroutine Particles_unitTest(fileUnit,perfect)

  use Particles_data, ONLY:  pt_numLocal, particles
  use Grid_interface, ONLY : Grid_getBlkBoundBox
  
#include "Flash.h"
#include "constants.h"

  implicit none

  
  integer, intent(IN)           :: fileUnit ! Output to file
  logical, intent(INOUT)        :: perfect  ! Flag to indicate errors

  integer :: p, blkID
  real,dimension(2,MDIM) :: limitBndBox

  

  do p = 1, pt_numLocal
     blkID=particles(BLK_PART_PROP,p)
     call Grid_getBlkBoundBox(blkID,limitBndBox)

     if(particles(POSX_PART_PROP,p)<limitBndBox(1,1))then
        write(fileUnit,*)'particle p = ',p,' out of limit lowerface iaxis', limitBndBox(1,1)
        perfect = .false.
     end if

     if(particles(POSX_PART_PROP,p)>limitBndBox(2,1)) then
        write(fileUnit,*)'particle p = ',p,' out of limit upperface iaxis',limitBndBox(2,1)
        perfect = .false.
     end if

     if(particles(POSY_PART_PROP,p)<limitBndBox(1,2)) then
        write(fileUnit,*)'particle p = ',p,' out of limit lowerface jaxis',limitBndBox(1,2)
        perfect = .false.  
     end if

     if(particles(POSY_PART_PROP,p)>limitBndBox(2,2)) then
        write(fileUnit,*)'particle p = ',p,' out of limit upperface jaxis',limitBndBox(2,2)
        perfect = .false.  
     end if

     if(particles(POSZ_PART_PROP,p)<limitBndBox(1,3)) then
        write(fileUnit,*)'particle p = ',p,' out of limit lowerface kaxis',limitBndBox(1,3)
        perfect = .false.     
     end if

     if(particles(POSZ_PART_PROP,p)>limitBndBox(2,3)) then
        write(fileUnit,*)'particle p = ',p,' out of limit upperface kaxis',limitBndBox(2,3)
        perfect = .false.    
     end if

  enddo
  
890 format(I5,5X,I5,5X,G12.4,3X,G15.7,3X,I5)
900 format(I5,3X,I6,2X,I5,3X,6(G12.5,3X))




  return

 
end subroutine Particles_unitTest
