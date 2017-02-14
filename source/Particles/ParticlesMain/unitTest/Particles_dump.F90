!!****if* source/Particles/ParticlesMain/unitTest/Particles_dump
!!
!! NAME
!!    Particles_dump
!!
!! SYNOPSIS
!!    Particles_dump()
!!                   integer(IN) :: blockCount,
!!                   integer(IN), dimension(:)  :: blockList,
!!                   integer(IN) :: nstep,
!!                   real(IN)    :: time)
!!
!! DESCRIPTION
!!
!!   Dump particle information to a plain old file.  The output file is called
!!      test_ParticlesDump_00##  where the ## relates to the processor number.
!!   Quick and dirty output implemented for blue gene testing.  You can see an
!!      example of the output by running setups/unitTest/Particles.
!!   However, there is an associated fidlr3/IDL routine particles_dump.pro which will
!!      read the input from this file.
!!
!! ARGUMENTS
!!
!!  
!!  blockCount:             integer(IN)     number of blocks within this processor
!!  blockList(blockCount):  integer(IN)     block IDs within this processor
!!  nstep:                  integer(IN)     current time step index
!!  time:                   real(IN)        current simulation time
!!
!! PARAMETERS
!!
!!***


!-------------------------------------------------------------------
subroutine Particles_dump(blockCount,blockList,nstep,time,dt)

  use Particles_data, ONLY:  pt_numLocal, particles, pt_meshMe
  use Grid_interface, ONLY : Grid_getBlkBoundBox

#include "Flash.h"
#include "constants.h"

  implicit none

  
  integer, intent(IN) :: blockCount
  integer, intent(IN) :: blockList(blockCount)
  integer, intent(IN) :: nstep
  real, intent(IN)    :: time, dt

  integer ::  fileUnit = 2
  character(len=23)   :: filename
  integer,dimension(4) :: prNum
  integer :: temp, i, p
  real,dimension(2,MDIM) :: limit

! some trick here to get a filename that reflects the processor or time step
  temp = pt_meshMe
  call Grid_getBlkBoundBox(1,limit)
  do i = 1,4
     prNum(i)= mod(temp,10)
     temp = temp/10
  end do
  filename = "test_ParticlesDump_"//char(48+prNum(4))//char(48+prNum(3))//&
                                 char(48+prNum(2))//char(48+prNum(1))

  if(nstep==0) then
     open(fileUnit,file=filename)
     write(fileUnit,'("P",I0)') pt_meshMe
  else
     open(fileUnit,file=filename,position='append')
  end if

! remember that you're writing out the velocity from the OLD time step
  write(fileUnit,890)pt_meshMe,nstep,time,dt,pt_numLocal
  do p = 1, pt_numLocal
     !for unit test only write out if particles fail

!!$     write(fileUnit,900)p, int(particles(TAG_PART_PROP,p)),int(particles(BLK_PART_PROP,p)), &
!!$          particles(POSX_PART_PROP,p),particles(POSY_PART_PROP,p),particles(POSZ_PART_PROP,p), &
!!$          particles(VELX_PART_PROP,p),particles(VELY_PART_PROP,p),particles(VELZ_PART_PROP,p)

     if(particles(POSX_PART_PROP,p)<limit(1,1))write(fileUnit,*)'particle p = ',p,' out of limit lowerface iaxis', limit(1,1)
     if(particles(POSX_PART_PROP,p)>limit(2,1))write(fileUnit,*)'particle p = ',p,' out of limit upperface iaxis',limit(2,1)
     if(particles(POSY_PART_PROP,p)<limit(1,2))write(fileUnit,*)'particle p = ',p,' out of limit lowerface jaxis',limit(1,2)
     if(particles(POSY_PART_PROP,p)>limit(2,2))write(fileUnit,*)'particle p = ',p,' out of limit upperface jaxis',limit(2,2)
     if(particles(POSZ_PART_PROP,p)<limit(1,3))write(fileUnit,*)'particle p = ',p,' out of limit lowerface kaxis',limit(1,3)
     if(particles(POSZ_PART_PROP,p)>limit(2,3))write(fileUnit,*)'particle p = ',p,' out of limit upperface kaxis',limit(2,3)


  enddo
  
  close(fileUnit)
      
890 format(I5,5X,I5,5X,G12.4,3X,G15.7,3X,I5)
900 format(I5,3X,I6,2X,I5,3X,6(G12.5,3X))
 
  
end subroutine Particles_dump

