!!****if* source/Particles/ParticlesInitialization/Particles_addNew
!!
!! NAME
!!    Particles_addNew
!!
!! SYNOPSIS
!!    call Particles_addNew( integer(in)  :: count,
!!                  optional,real(in)     :: pos(MDIM,count),
!!                           logical(out) :: success)
!!
!! DESCRIPTION
!!
!!    This routine allows particles to be added during evolution.
!!    In the particles data structure it always initializes the tag and
!!    processor ID. If the optional argument "pos" is present then it
!!    will also initialize the position and block ID attributes in the
!!    particles. It returns the value FALSE in success if there isn't
!!    enough space left in the data structure for the requested number
!!    of particles.
!!
!! ARGUMENTS
!!
!!     count   :: the count of particles to be added
!!     pos     :: optional, contains the coordinates of the particles
!!     success :: This arg returns TRUE if there was enough space 
!!                in the particles data structure to add the requested
!!                number of particles, FALSE otherwise.
!!
!!    
!!  NOTES
!!
!!   This routine must be called collectively, i.e., by all MPI tasks
!!   in the pt_meshComm communicator, if the optional "pos" argument
!!   is present.
!!
!!   The constant MDIM is defined in constants.h .
!!
!!***

!!#define DEBUG_PARTICLES

subroutine Particles_addNew (count, pos, success)
  
  use Particles_data, ONLY : particles, &
       pt_maxPerProc, pt_numLocal, pt_meshComm, pt_meshMe, pt_indexList, &
       pt_indexCount

  use Logfile_interface, ONLY : Logfile_stampMessage
  use Grid_interface, ONLY : Grid_moveParticles
  use pt_interface, ONLY : pt_findTagOffset
  implicit none
  
#include "constants.h"
#include "Flash.h"
#include "Particles.h"
#include "Flash_mpi.h"

  integer, INTENT(in) :: count
  real, optional, dimension(MDIM,count), intent(IN)::pos
  logical, intent(OUT) :: success

  integer :: i, tagOffset, ierr, effCount
  logical,parameter :: coords_in_blk=.true.
  logical :: doAll, doAny, doLocal(2),doGlobal(2)
  character(len=80) :: message
!-------------------------------------------------------------------------------
  doLocal(1) = ((pt_numLocal+count).le.pt_maxPerProc)
  doLocal(2) = (.NOT. doLocal(1))

  call MPI_AllReduce(doLocal(1),doGlobal(1),2,MPI_LOGICAL,MPI_LAND,pt_meshComm,ierr)
  doAll = doGlobal(1)
  doAny = .NOT. doGlobal(2)

!  if((pt_numLocal+count).le.pt_maxPerProc) then
  if(doAny) then
     if (doLocal(1)) then
        success=.true.
        effCount = count
     else
        success = (count==0)
        effCount = 0
     end if
     call pt_findTagOffset(effCount,tagOffset)
     do i = 1,effCount
        particles(PROC_PART_PROP,pt_numLocal+i)=pt_meshMe
        particles(TAG_PART_PROP,pt_numLocal+i)=tagOffset+i
     end do
     if(present(pos)) then
        do i = 1,effCount
           particles(POSX_PART_PROP:POSZ_PART_PROP,pt_numLocal+i)=&
                pos(IAXIS:KAXIS,i)
        end do
        
        pt_numLocal=pt_numLocal+effCount
        call Grid_moveParticles(particles,NPART_PROPS,pt_maxPerProc,&
             pt_numLocal,pt_indexList,pt_indexCount,coords_in_blk)
     else
        pt_numLocal=pt_numLocal+effCount
     end if
     if (.NOT. success) then
98      format('Particles_addNew for',I6,' particles failed')
        write(message,98) count
        print *,trim(message),' on',pt_meshMe
        call Logfile_stampMessage(message, force=.TRUE.)
     end if
  else
99   format('Particles_addNew for',I6,' particles failed everywhere') 
     write(message,99) count
     print *,trim(message),', detected on',pt_meshMe
     call Logfile_stampMessage(message, force=.TRUE.)
     success=.false.
  end if
  
  return
  
end subroutine Particles_addNew
