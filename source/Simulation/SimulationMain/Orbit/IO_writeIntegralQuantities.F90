!!****if* source/Simulation/SimulationMain/Orbit/IO_writeIntegralQuantities
!!
!!  NAME
!!    IO_writeIntegralQuantities
!!
!!  SYNOPSIS
!!    call IO_writeIntegralQuantities()
!!                                    integer(in) :: isFirst,
!!                                    real(in)    :: simTime)
!!
!!  DESCRIPTION
!!     Compute the values of integral quantities (eg. total energy)
!!     and write them to an ASCII file.  If this is the initial step,
!!     create the file and write a header to it before writing the
!!     data.
!!
!!     Presently, this supports 1, 2, and 3-d Cartesian geometry and 2-d
!!     cylindrical geometry (r,z).  More geometries can be added by modifying
!!     the volume of each zone (dvol) computed in SumLocalIntegrals.
!!
!!     Special version for particle orbit problem:  for each particle,
!!     writes out particle position, velocitie, and acceleration.  Each
!!     particle gets a separate file named 'pdata.0000', where 0000 is
!!     replaced by the particle's tag number.
!!
!!  ARGUMENTS
!!    
!!   isFirst - if 1 then write header info plus data, otherwise just write data
!!   simTime - simulation time
!!
!!
!!***

subroutine IO_writeIntegralQuantities ( isfirst, simtime)



  use Particles_data, ONLY : pt_numLocal, particles

   use IO_data, ONLY : io_globalMe
  implicit none

#include "Flash.h"
#include "Flash_mpi.h"
#include "constants.h"

  integer, intent(in) :: isfirst
  real, intent(in)              :: simtime

  character(len=80) :: outfile = "pdata.0000"
  integer           :: i, ierr
  logical, save     :: firstCall = .TRUE.

  real :: K, W, M, Px, Py, Pz
  real :: myK, myW, myM, myPx, myPy, myPz
  real :: vx, vy, vz, mass, gpot

!=============================================================================

  if (firstCall) then
     if (pt_numLocal > 0) then
        do i = 1, pt_numLocal
           write (outfile(7:10),'(I4.4)') int(particles(TAG_PART_PROP,i))
           open (11, file=outfile, status="replace")
           close (11)
        enddo
     endif

     !Request that only the master processor creates flash.dat.
     if (io_globalMe  == 0) then
        open (11, file='flash.dat', status="replace")
        close(11)
     end if

     firstCall = .FALSE.
  endif


! Write positions, velocities, and accelerations for active particles

  myK  = 0.
  myW  = 0.
  myPx = 0.
  myPy = 0.
  myPz = 0.
  myM  = 0.

  if (pt_numLocal > 0) then
     do i = 1, pt_numLocal
        write (outfile(7:10),'(I4.4)') int(particles(TAG_PART_PROP,i))
        open (11, file=outfile, position="append")
        write (11,999) simtime, particles(POSX_PART_PROP,i), &
             particles(POSY_PART_PROP,i), &
             particles(POSZ_PART_PROP,i), &
             particles(VELX_PART_PROP,i), &
             particles(VELY_PART_PROP,i), &
             particles(VELZ_PART_PROP,i), &
             particles(ACCX_PART_PROP,i), &
             particles(ACCY_PART_PROP,i), &
             particles(ACCZ_PART_PROP,i)
        close (11)
        vx = particles(VELX_PART_PROP,i)
        vy = particles(VELY_PART_PROP,i)
        vz = particles(VELZ_PART_PROP,i)
        mass = particles(MASS_PART_PROP,i)
        gpot = particles(GPOT_PART_PROP,i)
        myK  = myK  + 0.5*mass*(vx**2+vy**2+vz**2)
        myW  = myW  + 0.5*mass*gpot
        myPx = myPx + mass*vx
        myPy = myPy + mass*vy
        myPz = myPz + mass*vz
        myM  = myM  + mass
     enddo
  endif

  call mpi_allreduce(myM,  M,  1, FLASH_REAL, MPI_SUM, MPI_COMM_WORLD, ierr)
  call mpi_allreduce(myK,  K,  1, FLASH_REAL, MPI_SUM, MPI_COMM_WORLD, ierr)
  call mpi_allreduce(myW,  W,  1, FLASH_REAL, MPI_SUM, MPI_COMM_WORLD, ierr)
  call mpi_allreduce(myPx, Px, 1, FLASH_REAL, MPI_SUM, MPI_COMM_WORLD, ierr)
  call mpi_allreduce(myPy, Py, 1, FLASH_REAL, MPI_SUM, MPI_COMM_WORLD, ierr)
  call mpi_allreduce(myPz, Pz, 1, FLASH_REAL, MPI_SUM, MPI_COMM_WORLD, ierr)
  
 !Request that only the master processor writes to flash.dat.
  if (io_globalMe  == 0) then
     open(11, file='flash.dat', position="append")
     write(11,998) simtime, M, K, W, K+W, Px, Py, Pz
     close(11)
  endif

998 format (8(ES12.5, :, 1X))
999 format (10(ES12.5, :, 1X))

!=============================================================================

  return
end subroutine IO_writeIntegralQuantities
