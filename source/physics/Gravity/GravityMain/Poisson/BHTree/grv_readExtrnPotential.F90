!!****if* source/physics/Gravity/GravityMain/Poisson/BHTree/grv_readExtrnPotential
!!
!! NAME
!!
!!  grv_readExtrnPotential
!!
!!
!! SYNOPSIS
!!
!!  call grv_readExtrnPotential()
!!
!! DESCRIPTION
!!
!!   Reads external gravitational field from file (filename set by runtime 
!!   parameter grv_bhExtrnPotFile) and stores it into arrays
!!   grv_bhExtrnPotCoord, grv_bhExtrnPotPot, grv_bhExtrnPotAcc.
!!
!! ARGUMENTS
!!
!!
!! RESULT
!!
!!
!!***

subroutine grv_readExtrnPotential()
  use Gravity_data, ONLY : grv_meshMe, grv_meshComm, grv_bhExtrnPotNMax, &
    grv_bhExtrnPotCoord, grv_bhExtrnPotPot, grv_bhExtrnPotAcc, &
    grv_bhExtrnPotDel, grv_useExternalPotential, grv_bhExtrnPotIType, &
    grv_bhEPTypeR, grv_bhEPTypeX, grv_bhEPTypeY, grv_bhEPTypeZ, &
    grv_bhExtrnPotType, grv_bhExtrnPotFile, grv_bhExtrnPotN
  use Driver_interface, ONLY : Driver_abortFlash
  implicit none
#include "constants.h"
#include "Flash_mpi.h"
  integer :: i, istat, io_status
  real :: readArr(3,grv_bhExtrnPotNMax)

  if (.not. grv_useExternalPotential) return

  if (grv_meshMe .eq. MASTER_PE) then
    open(unit = 55, file = grv_bhExtrnPotFile, status = 'old', iostat = io_status)
    if (io_status /= 0) then
      print *, "Unable to open External potential file: ", grv_bhExtrnPotFile
      call Driver_abortFlash("Could not open external potential file.")
    endif
    i = 1
    do
      read(55,*, iostat = io_status) readArr(:,i)
      if (io_status /= 0) exit
      i = i + 1
      if (i > grv_bhExtrnPotNMax) &
        & call Driver_abortFlash("Too many lines in the external potential file.")
    enddo
    close(unit = 55)
    grv_bhExtrnPotN = i - 1
  endif
  call MPI_Bcast(grv_bhExtrnPotN, 1, FLASH_INTEGER, MASTER_PE, grv_meshComm, istat)

  ! then, read and communicate radial coordinate and radial grav force
  allocate(grv_bhExtrnPotCoord(grv_bhExtrnPotN), stat=istat)
  if (istat .ne. 0) call Driver_abortFlash("Could not allocate grv_bhExtrnPotCoord")
  allocate(grv_bhExtrnPotPot(grv_bhExtrnPotN), stat=istat)
  if (istat .ne. 0) call Driver_abortFlash("Could not allocate grv_bhExtrnPotPot")
  allocate(grv_bhExtrnPotAcc(grv_bhExtrnPotN), stat=istat)
  if (istat .ne. 0) call Driver_abortFlash("Could not allocate grv_bhExtrnPotAcc")
  grv_bhExtrnPotCoord = readArr(1,1:grv_bhExtrnPotN)
  grv_bhExtrnPotPot   = readArr(2,1:grv_bhExtrnPotN)
  grv_bhExtrnPotAcc   = readArr(3,1:grv_bhExtrnPotN)
  call MPI_Bcast(grv_bhExtrnPotCoord, grv_bhExtrnPotN, FLASH_REAL, MASTER_PE, grv_meshComm, istat)
  call MPI_Bcast(grv_bhExtrnPotPot,   grv_bhExtrnPotN, FLASH_REAL, MASTER_PE, grv_meshComm, istat)
  call MPI_Bcast(grv_bhExtrnPotAcc,   grv_bhExtrnPotN, FLASH_REAL, MASTER_PE, grv_meshComm, istat)

  grv_bhExtrnPotDel = (grv_bhExtrnPotCoord(grv_bhExtrnPotN) - grv_bhExtrnPotCoord(1)) &
  &                 / (grv_bhExtrnPotN - 1)

  if (grv_meshMe .eq. MASTER_PE) then
    print *, "Gravity external potential read."
    print *, "  grv_bhExtrnPotN = ", grv_bhExtrnPotN
    print *, "  grv_bhExtrnPotDel = ", grv_bhExtrnPotDel
    !do i = 1, grv_bhExtrnPotN
    !  print *, "    i, coord, pot, acc = ", i, grv_bhExtrnPotCoord(i) &
    !  , grv_bhExtrnPotPot(i), grv_bhExtrnPotAcc(i)
    !enddo
  endif


  select case (grv_bhExtrnPotType)
    case("spherical")
      grv_bhExtrnPotIType = grv_bhEPTypeR
    case("planex")
      grv_bhExtrnPotIType = grv_bhEPTypeX
      call Driver_abortFlash("grv_readExtrnPotential: planex symmetry not supported yet")
    case("planey")
      grv_bhExtrnPotIType = grv_bhEPTypeY
      call Driver_abortFlash("grv_readExtrnPotential: planey symmetry not supported yet")
    case("planez")
      grv_bhExtrnPotIType = grv_bhEPTypeZ
    case default
      call Driver_abortFlash("grv_readExtrnPotential: unrecognized potential type")
  end select

end subroutine grv_readExtrnPotential


