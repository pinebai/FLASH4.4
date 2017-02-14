!!****if* source/Simulation/SimulationMain/RTFlame/IO_writeIntegralQuantities
!!
!! NAME
!!
!!  IO_writeIntegralQuantities
!!
!! SYNOPSIS
!!
!!  call IO_writeIntegralQuantities(integer(in) :: isfirst,
!!                                  real(in) :: simtime)
!!
!! DESCRIPTION
!!
!!
!! ARGUMENTS
!!
!!   isfirst : 
!!
!!   simtime : 
!!
!! AUTOGENROBODOC
!!
!!
!!***

! See source/IO/IOMain/IO_writeIntegralQuantities.F90
! for API and original (example) subroutine
!
! This version has some additional metrics specific to the RT Flame
! in a channel (burning rate, surface area)
!
! Dean Townsley 2008


!!REORDER(4):solnData

subroutine IO_writeIntegralQuantities (isFirst, simTime)

  use IO_data, ONLY : io_restart, io_statsFileName, io_globalMe, io_globalComm
  use Grid_interface, ONLY : Grid_getListOfBlocks, &
    Grid_getBlkIndexLimits, Grid_getBlkPtr, Grid_getSingleCellVol, &
    Grid_releaseBlkPtr, Grid_fillGuardCells
  use IO_interface, ONLY : IO_getScalar, IO_setScalar
  use Driver_data, ONLY : dr_dtOld
  use Simulation_data, ONLY : sim_dens_u, sim_last_burned_mass
  use ut_contourSurfaceInterface, ONLY: ut_contourSurfaceAreaBlock
  use fl_fsData, ONLY : fl_fsConstFlameSpeed

  implicit none

#include "Flash_mpi.h"
#include "constants.h"
#include "Flash.h"
  
  real, intent(in) :: simTime

  integer, intent(in) :: isFirst

  integer :: lb, count
  
  integer :: funit = 99
  integer :: error
  
  character (len=MAX_STRING_LENGTH), save :: fname 
  
  integer :: blockList(MAXBLOCKS)

  integer :: blkLimits(HIGH, MDIM), blkLimitsGC(HIGH, MDIM)

  integer, parameter ::  nGlobalSum = 14  ! Number of globally-summed quantities
  real :: gsum(nGlobalSum) !Global summed quantities
  real :: lsum(nGlobalSum) !Global summed quantities

  integer :: i, j, k
  real :: dvol             !, del(MDIM)
  real, DIMENSION(:,:,:,:), POINTER :: solnData

  integer,parameter :: nlevels = 3
  real, dimension(nlevels) :: isolevels, blkAreas
  integer :: point(MDIM)
  real    :: dt, brate, fspd_ratio, fspd

  logical :: gcMask(NUNK_VARS)

  isolevels(1) = 0.1
  isolevels(2) = 0.5
  isolevels(3) = 0.9

#if NDIM > 1
  gcMask(:) = .false.
  gcMask(FLAM_MSCALAR) = .true.

  call Grid_fillGuardCells(CENTER, ALLDIR,&
       maskSize=NUNK_VARS, mask=gcMask,makeMaskConsistent=.true.)
#endif
  

  ! Sum quantities over all locally held leaf-node blocks.
  gsum  = 0.
  lsum = 0.
  
  call Grid_getListOfBlocks(LEAF, blockList, count)
  
  do lb = 1, count
     !get the index limits of the block
     call Grid_getBlkIndexLimits(blockList(lb), blkLimits, blkLimitsGC)

     ! get a pointer to the current block of data
     call Grid_getBlkPtr(blockList(lb), solnData)

     ! Sum contributions from the indicated blkLimits of cells.
     do k = blkLimits(LOW,KAXIS), blkLimits(HIGH,KAXIS)
        do j = blkLimits(LOW,JAXIS), blkLimits(HIGH,JAXIS)
           do i = blkLimits(LOW,IAXIS), blkLimits(HIGH,IAXIS)
              
              point(IAXIS) = i
              point(JAXIS) = j
              point(KAXIS) = k

!! Get the cell volume for a single cell
              call Grid_getSingleCellVol(blockList(lb), EXTERIOR, point, dvol)
     
              ! mass   
#ifdef DENS_VAR
              lsum(1) = lsum(1) + solnData(DENS_VAR,i,j,k)*dvol 
#endif           


#ifdef DENS_VAR
#ifdef VELX_VAR      
              ! momentum
              lsum(2) = lsum(2) + solnData(DENS_VAR,i,j,k) * & 
                   &                                solnData(VELX_VAR,i,j,k)*dvol 
           
#endif
#ifdef VELY_VAR      

              lsum(3) = lsum(3) + solnData(DENS_VAR,i,j,k) * & 
                   &                                solnData(VELY_VAR,i,j,k)*dvol
           
#endif
#ifdef VELZ_VAR      
              lsum(4) = lsum(4) + solnData(DENS_VAR,i,j,k) * & 
                   &                                solnData(VELZ_VAR,i,j,k)*dvol
#endif

              ! total energy
#ifdef ENER_VAR
              lsum(5) = lsum(5) + solnData(ENER_VAR,i,j,k) * & 
                   &                                solnData(DENS_VAR,i,j,k)*dvol
#endif
           
#ifdef VELX_VAR      
#ifdef VELY_VAR      
#ifdef VELZ_VAR      
              ! kinetic energy
              lsum(6) = lsum(6) + 0.5*solnData(DENS_VAR,i,j,k) * & 
                   &                             (solnData(VELX_VAR,i,j,k)**2+ & 
                   &                              solnData(VELY_VAR,i,j,k)**2+ & 
                   &                              solnData(VELZ_VAR,i,j,k)**2)*dvol           

#endif
#endif
#endif

#ifdef TURB_VAR
              ! turbulent kinetic energy
              lsum(7) = lsum(7) + 0.5*solnData(DENS_VAR,i,j,k) * &
                   &                  solnData(TURB_VAR,i,j,k)**2 * dvol
#endif

#ifdef EINT_VAR
              ! internal energy
              lsum(8) = lsum(8) + solnData(DENS_VAR,i,j,k) * & 
                   &                                solnData(EINT_VAR,i,j,k)*dvol
#endif
#endif ! ifdef DENS_VAR

              lsum(9) = lsum(9) + solnData(FLAM_MSCALAR,i,j,k) * &
                                         solnData(DENS_VAR,i,j,k)*dvol

              if (solnData(FLAM_MSCALAR,i,j,k) > 1.e-6 .and.  solnData(FLAM_MSCALAR,i,j,k) < 1.e-3) then
                 ! find average density just ahead of flame
                 lsum(10) = lsum(10) + solnData(DENS_VAR,i,j,k)*dvol
                 lsum(11) = lsum(11) + dvol
              endif

           enddo
        enddo
     enddo

#if NDIM > 1
     call ut_contourSurfaceAreaBlock(nlevels,isolevels,solnData(FLAM_MSCALAR,:,:,:), &
                                     blkLimits,blockList(lb),blkAreas)
     lsum(12) = lsum(12) + blkAreas(1)
     lsum(13) = lsum(13) + blkAreas(2)
     lsum(14) = lsum(14) + blkAreas(3)
#endif


     call Grid_releaseBlkPtr(blockList(lb), solnData)

  enddo
  

  
  ! Now the MASTER_PE sums the local contributions from all of
  ! the processors and writes the total to a file.
  
  call MPI_Reduce (lsum, gsum, nGlobalSum, FLASH_REAL, MPI_SUM, & 
       &                MASTER_PE, io_globalComm, error)
  
  dt = dr_dtOld
  if (io_globalMe == MASTER_PE) then
     ! gsum only valid on master node
     if (sim_last_burned_mass==-1.0) then
           brate=0.0
     else
        brate = (gsum(9)-sim_last_burned_mass)/(2*dt)
     endif
     sim_last_burned_mass = gsum(9)
  else
     sim_last_burned_mass = 0.0
  endif
  ! has to go in the linked list on all nodes for IO to work right
  call IO_setScalar("last_burned_mass",sim_last_burned_mass)

  if (io_globalMe == MASTER_PE) then
     ! average density just ahead of flame
     if (gsum(11) > 0.0) then
        gsum(10) = gsum(10)/gsum(11)
     else
        gsum(10) = sim_dens_u
     endif

     ! calculate flame front propagation speed: s = mb / rho / BurnArea
     fspd = brate / gsum(10) / gsum(12)
     fspd_ratio = fspd / fl_fsConstFlameSpeed
     
     ! create the file from scratch if it is a not a restart simulation, 
     ! otherwise append to the end of the file
     if (isfirst == 0) then
        open (funit, file=trim(io_statsFileName), position='APPEND')
     else 
        if (.NOT. io_restart) then
           open (funit, file=trim(io_statsFileName)) 
           write (funit, 10)               &
                '#time                     ', &
                'mass                      ', &
                'x-momentum                ', &
                'y-momentum                ', & 
                'z-momentum                ', &
                'E_total                   ', &
                'E_kinetic                 ', &
                'E_turbulent               ', &
                'E_internal                ', &
                'Burned Mass               ', &
                'dens_burning_ave          ', &
                'db_ave samplevol          ', &
                'Burning rate              ', &
                'fspd to input_fspd ratio  ', &
                'surface area flam=0.1     ', &
                'surface area flam=0.5     ', &
                'surface area flam=0.9     '

10         format (2x,50(a25, :, 1X))

        else
           open (funit, file=trim(io_statsFileName), position='APPEND')
           write (funit, 11) 
11         format('# simulation restarted')
        endif
     endif

     ! Write the global sums to the file.
     write (funit, 12) simtime, gsum(1:11),brate,fspd_ratio,gsum(12:14)
12   format (1x, 50(es25.18, :, 1x))
 
     close (funit)          ! Close the file.
     
  endif
  
  call MPI_Barrier (io_globalComm, error)
  
  !=============================================================================
  
  return
end subroutine IO_writeIntegralQuantities



