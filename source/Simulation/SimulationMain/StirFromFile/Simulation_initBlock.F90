!!****if* source/Simulation/SimulationMain/StirFromFile/Simulation_initBlock
!!
!! NAME
!!  Simulation_initBlock
!!
!! SYNOPSIS
!!  Simulation_initBlock(integer, intent(IN)  :: blockid)
!!
!! DESCRIPTION
!!  Initializes data (density, pressure, velocity, etc.) for
!!  a specified block. 
!!
!! ARGUMENTS
!!   blockid : ID of block in current processor
!!
!! AUTHOR
!!   Christoph Federrath, 2008-2013
!!
!!***


subroutine Simulation_initBlock (blockID)

  use Driver_interface, ONLY : Driver_getSimTime, Driver_getMype
  use Simulation_data, ONLY  : sim_rhoAmbient, sim_cAmbient, sim_gamma, &
                               sim_magnetic, sim_MagField_z
  use Grid_interface, ONLY   : Grid_getBlkIndexLimits, Grid_getCellCoords, &
                               Grid_getBlkPtr, Grid_releaseBlkPtr
  use Eos_interface, ONLY    : Eos_wrapped

  implicit none 

#include "constants.h"
#include "Flash.h"

  integer, intent(IN) :: blockID

  integer, save                :: myPE
  integer                      :: sizeZ, sizeY, sizeX, i, j, k
  integer, dimension(2,MDIM)   :: blkLimits, blkLimitsGC
  logical, parameter           :: gcell = .true.
  real                         :: del(MDIM), dvol

  real, DIMENSION(:,:,:,:), POINTER :: solnData
#if NFACE_VARS > 0
  real, DIMENSION(:,:,:,:), POINTER :: facexData, faceyData, facezData
#endif
#ifdef FIXEDBLOCKSIZE
  real, dimension(GRID_IHI_GC) :: xCoord
  real, dimension(GRID_JHI_GC) :: yCoord
  real, dimension(GRID_KHI_GC) :: zCoord
#else
  real, allocatable, dimension(:) :: xCoord, yCoord, zCoord
  integer :: istat
#endif

  logical, parameter :: Debug = .false.


  call Driver_getMype(GLOBAL_COMM, myPE)

  if (Debug) print *, '[', myPE, '] Simulation_initBlock entering.'

  !get the index limits of the block
  call Grid_getBlkIndexLimits(blockID, blkLimits, blkLimitsGC)

  ! get a pointer to the current block of data
  call Grid_getBlkPtr(blockID, solnData)
  if (Debug) print *, '[', myPE, '] Simulation_initBlock: got solnData pointer.'

  !getting the dx's
  call Grid_getDeltas(blockID, del)

#if NDIM == 1
   dvol = del(IAXIS)
#endif
#if NDIM == 2
   dvol = del(IAXIS) * del(JAXIS)
#endif
#if NDIM == 3
   dvol = del(IAXIS) * del(JAXIS) * del(KAXIS)
#endif

   sizeX = blkLimitsGC(HIGH,IAXIS)-blkLimitsGC(LOW,IAXIS)+1
   sizeY = blkLimitsGC(HIGH,JAXIS)-blkLimitsGC(LOW,JAXIS)+1
   sizeZ = blkLimitsGC(HIGH,KAXIS)-blkLimitsGC(LOW,KAXIS)+1

#ifndef FIXEDBLOCKSIZE
   allocate(xCoord(sizeX),stat=istat)
   if (istat .ne. 0) call Driver_abortFlash("could not allocate xCoord in Simulation_initBlock.F90")
   allocate(yCoord(sizeY),stat=istat)
   if (istat .ne. 0) call Driver_abortFlash("could not allocate yCoord in Simulation_initBlock.F90")
   allocate(zCoord(sizeZ),stat=istat)
   if (istat .ne. 0) call Driver_abortFlash("could not allocate zCoord in Simulation_initBlock.F90")
#endif

   ! x coordinates
   call Grid_getCellCoords(IAXIS,blockID,CENTER,gcell,xCoord,sizeX)
#if NDIM > 1
   ! y coordinates
   call Grid_getCellCoords(JAXIS,blockID,CENTER,gcell,yCoord,sizeY)
#if NDIM > 2
   ! z coordinates
   call Grid_getCellCoords(KAXIS,blockID,CENTER,gcell,zCoord,sizeZ)
#endif
#endif

   ! loop over cells in block
   do k = blkLimits(LOW,KAXIS), blkLimits(HIGH,KAXIS)
       do j = blkLimits(LOW,JAXIS), blkLimits(HIGH,JAXIS)
           do i = blkLimits(LOW,IAXIS), blkLimits(HIGH,IAXIS)

              solnData(DENS_VAR,i,j,k) = sim_rhoAmbient
              solnData(VELX_VAR,i,j,k) = 0.
              solnData(VELY_VAR,i,j,k) = 0.
              solnData(VELZ_VAR,i,j,k) = 0.
              solnData(PRES_VAR,i,j,k) = solnData(DENS_VAR,i,j,k)*sim_cAmbient**2
              solnData(EINT_VAR,i,j,k) = solnData(PRES_VAR,i,j,k)/((sim_gamma-1.0)*solnData(DENS_VAR,i,j,k))
              solnData(ENER_VAR,i,j,k) = solnData(EINT_VAR,i,j,k) + &
                                         0.5*( solnData(VELX_VAR,i,j,k)**2 + &
                                               solnData(VELY_VAR,i,j,k)**2 + &
                                               solnData(VELZ_VAR,i,j,k)**2   )

              if (sim_magnetic) then
#ifdef MAGX_VAR
                 solnData(MAGX_VAR,i,j,k) = 0.0
#endif
#ifdef MAGY_VAR
                 solnData(MAGY_VAR,i,j,k) = 0.0
#endif
#ifdef MAGZ_VAR
                 solnData(MAGZ_VAR,i,j,k) = sim_MagField_z
#endif
#ifdef MAGP_VAR
                 solnData(MAGP_VAR,i,j,k) = 0.5*( solnData(MAGX_VAR,i,j,k)**2 + &
                                                  solnData(MAGY_VAR,i,j,k)**2 + &
                                                  solnData(MAGZ_VAR,i,j,k)**2   )
#endif
              endif

           enddo ! i
       enddo ! j
   enddo ! k

   call Grid_releaseBlkPtr(blockID, solnData)

#ifndef FIXEDBLOCKSIZE
   deallocate(xCoord)
   deallocate(yCoord)
   deallocate(zCoord)
#endif

#if NFACE_VARS > 0
   ! in case we are running the unsplit staggered mesh solver, initialize the B field also at faces
   if (sim_magnetic) then

    if (Debug) print *, '[', myPE, '] --- initializing face variables for magnetic field components ...'
    if (Debug) print *, '[', myPE, '] getting pointer to solndata ...'
    call Grid_getBlkPtr(blockID,solnData,CENTER)
    if (Debug) print *, '[', myPE, '] getting pointer to facexData ...'
    call Grid_getBlkPtr(blockID,facexData,FACEX)
    if (Debug) print *, '[', myPE, '] getting pointer to faceyData ...'
    call Grid_getBlkPtr(blockID,faceyData,FACEY)
    if (Debug) print *, '[', myPE, '] getting pointer to facezData ...'
    call Grid_getBlkPtr(blockID,facezData,FACEZ)
    if (Debug) print *, '[', myPE, '] ... getting pointers done.'

    ! Loop over cells in the block for face magnetic field values
    do k = blkLimitsGC(LOW,KAXIS),blkLimitsGC(HIGH,KAXIS)
       do j = blkLimitsGC(LOW,JAXIS),blkLimitsGC(HIGH,JAXIS)
          do i = blkLimitsGC(LOW,IAXIS),blkLimitsGC(HIGH,IAXIS)
              ! Cell face-centered variables for StaggeredMesh scheme
              facexData(MAG_FACE_VAR,i,j,k) = 0.0
              faceyData(MAG_FACE_VAR,i,j,k) = 0.0
              facezData(MAG_FACE_VAR,i,j,k) = sim_MagField_z
          enddo
       enddo
    enddo

    i = blkLimitsGC(HIGH,IAXIS) + 1
    do k = blkLimitsGC(LOW,KAXIS), blkLimitsGC(HIGH,KAXIS)
      do j = blkLimitsGC(LOW,JAXIS), blkLimitsGC(HIGH,JAXIS)
        facexData(MAG_FACE_VAR,i,j,k) = 0.0
      enddo
    enddo

    j = blkLimitsGC(HIGH,JAXIS) + 1
    do k = blkLimitsGC(LOW,KAXIS), blkLimitsGC(HIGH,KAXIS)
      do i = blkLimitsGC(LOW,IAXIS), blkLimitsGC(HIGH,IAXIS)
        faceyData(MAG_FACE_VAR,i,j,k) = 0.0
      enddo
    enddo

    k = blkLimitsGC(HIGH,KAXIS) + 1
    do j = blkLimitsGC(LOW,JAXIS), blkLimitsGC(HIGH,JAXIS)
      do i = blkLimitsGC(LOW,IAXIS), blkLimitsGC(HIGH,IAXIS)
        facezData(MAG_FACE_VAR,i,j,k) = sim_MagField_z
      enddo
    enddo

    ! Release pointer
    call Grid_releaseBlkPtr(blockID,facexData,FACEX)
    call Grid_releaseBlkPtr(blockID,faceyData,FACEY)
    call Grid_releaseBlkPtr(blockID,facezData,FACEZ)
    call Grid_releaseBlkPtr(blockID,solnData,CENTER)

    if (Debug) print *, '[', myPE, '] --- initialization finished.'

   endif
#endif
! end USM magnetic

  ! update eint and temp by calling EOS
  call Eos_wrapped(MODE_DENS_PRES, blkLimits, blockID)

  if (Debug) print *, '[', myPE, '] Simulation_initBlock exiting.'

end subroutine Simulation_initBlock
