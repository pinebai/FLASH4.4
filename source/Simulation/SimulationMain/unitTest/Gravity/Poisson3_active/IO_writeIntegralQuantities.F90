!!****if* source/Simulation/SimulationMain/unitTest/Gravity/Poisson3_active/IO_writeIntegralQuantities
!!
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
!!
!!   Compute the values of integral quantities (eg. total energy)
!!   and write them to an ASCII file.  If this is the initial step,
!!   create the file and write a header to it before writing the data.
!!
!!   Presently, this supports 1, 2, and 3-d Cartesian geometry and 2-d
!!   cylindrical geometry (r,z).  More geometries can be added by
!!   modifying the volume of each zone (dvol).
!!
!!   Users should modify this routine if they want to store any
!!   quantities other than default values in the flash.dat.  Make sure
!!   to modify the nGlobalSum parameter to match the number of
!!   quantities written.  Also make sure to modify the header to match
!!   the names of quantities with those calculated in the lsum and
!!   gsum arrays.
!!  
!!  ARGUMENTS
!!    
!!   isFirst - if 1 then write header info plus data, otherwise just write data
!!   simTime - simulation time
!!
!!
!!***

!!REORDER(4):solnData

subroutine IO_writeIntegralQuantities ( isFirst, simTime)

  use IO_data, ONLY : io_restart, io_statsFileName
  use Grid_interface, ONLY : Grid_getListOfBlocks, &
    Grid_getBlkIndexLimits, Grid_getBlkPtr, Grid_getSingleCellVol, &
    Grid_releaseBlkPtr, Grid_getCellCoords
  use Simulation_data, ONLY: sim_xctr, sim_yctr, sim_zctr, sim_a1, sim_a3, &
       sim_a1inv, sim_a3inv, sim_initGeometry

   use IO_data, ONLY : io_globalMe
  implicit none

#include "mpif.h"
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
  real, dimension(:), allocatable :: xCenter, yCenter, zCenter
  integer     :: xSize, ySize, zSize
  real        :: radius2, xdist, ydist, zdist
  real        :: radiusInside2, a3x90, a3x90Inv

  integer, parameter ::  nGlobalSum = 11          ! Number of globally-summed quantities
  real :: gsum(nGlobalSum) !Global summed quantities
  real :: lsum(nGlobalSum) !Global summed quantities

  integer :: i, j, k, ii, jj, kk
  real :: dvol             !, del(MDIM)
  real, DIMENSION(:,:,:,:), POINTER :: solnData

  integer :: point(MDIM)

  ! define a variable that is 90% of the polar radius
  a3x90 = 0.9*sim_a3
  a3x90Inv = 0.9*sim_a3inv

  ! Sum quantities over all locally held leaf-node blocks.
  gsum  = 0.
  lsum = 0.
  
  call Grid_getListOfBlocks(LEAF, blockList, count)
  
  do lb = 1, count
     !get the index limits of the block
     call Grid_getBlkIndexLimits(blockList(lb), blkLimits, blkLimitsGC)

     xSize = blkLimits(HIGH,IAXIS) - blkLimits(LOW,IAXIS) + 1
     ySize = blkLimits(HIGH,JAXIS) - blkLimits(LOW,JAXIS) + 1
     zSize = blkLimits(HIGH,KAXIS) - blkLimits(LOW,KAXIS) + 1
     allocate(xCenter(xSize))
     allocate(yCenter(ySize))
     allocate(zCenter(zSize))

     ! get a pointer to the current block of data
     call Grid_getBlkPtr(blockList(lb), solnData)

     !! Get the cell coordinates
     call Grid_getCellCoords(IAXIS,blockList(lb),CENTER, .false., xCenter, xSize)
     call Grid_getCellCoords(JAXIS,blockList(lb),CENTER, .false., yCenter, ySize)
     call Grid_getCellCoords(KAXIS,blockList(lb),CENTER, .false., zCenter, zSize)


     ! Sum contributions from the indicated blkLimits of cells.
     do k = blkLimits(LOW,KAXIS), blkLimits(HIGH,KAXIS)
        kk = k - blkLimits(LOW,KAXIS) + 1
        zdist = (zCenter(kk) - sim_zctr)
        do j = blkLimits(LOW,JAXIS), blkLimits(HIGH,JAXIS)
           jj = j - blklimits(LOW,JAXIS) + 1
           ydist = (yCenter(jj) - sim_yctr)
           do i = blkLimits(LOW,IAXIS), blkLimits(HIGH,IAXIS)
              ii = i - blkLimits(LOW,IAXIS) + 1
              xdist = (xCenter(ii) - sim_xctr)


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


#ifdef EINT_VAR
              ! internal energy
              lsum(7) = lsum(7) + solnData(DENS_VAR,i,j,k) * & 
                   &                                solnData(EINT_VAR,i,j,k)*dvol
#endif
#endif ! ifdef DENS_VAR

!! Specialist information for the Maclaurin spheriod case

              lsum(8) = lsum(8) + solnData(ERRN_VAR,i,j,k)**2

              !! Divide inside and outside spheriod, as well as inside a sphere 
              !!   of radius 0.9*a3 (completely inside sphereoid)
              select case (sim_initGeometry)
              case (CYLINDRICAL)   ! 2d axisymmetric
                 radius2 = (xdist*sim_a1inv)**2 + (ydist*sim_a3inv)**2
                 radiusInside2 = (xdist*a3x90Inv)**2 + (ydist*a3x90Inv)**2
              case (CARTESIAN)       ! 3d cartesian
                 radius2 = (xdist*sim_a1inv)**2 + (ydist*sim_a1inv)**2 + (zdist*sim_a3inv)**2
                 radiusInside2 = (xdist*a3x90Inv)**2 + (ydist*a3x90Inv)**2 + (zdist*a3x90Inv)**2
              end select
              if (radius2 < 1.0) then     !! inside the spheroid
                 lsum(9) = lsum(9) + solnData(ERRN_VAR,i,j,k)**2
              else                        !! outside the spheriod
                 lsum(10) = lsum(10) + solnData(ERRN_VAR,i,j,k)**2
              end if
              if (radiusInside2 < 1.0) then   !! completely enclosed by spheriod
                 lsum(11) = lsum(11) + solnData(ERRN_VAR,i,j,k)**2
              end if
           enddo
        enddo
     enddo

     call Grid_releaseBlkPtr(blockList(lb), solnData)
     deallocate(xCenter,yCenter,zCenter)

  enddo
  
 
  ! Now the MASTER_PE sums the local contributions from all of
  ! the processors and writes the total to a file.
  
  call MPI_Reduce (lsum, gsum, nGlobalSum, MPI_Double_Precision, MPI_Sum, & 
       &                MASTER_PE, MPI_Comm_World, error)
  
  ! Take sqrts for L2 norms
  lsum(8) = sqrt(lsum(8))
  lsum(9) = sqrt(lsum(9))
  lsum(10) = sqrt(lsum(10))
  lsum(11) = sqrt(lsum(11))

  if (io_globalMe  == MASTER_PE) then
     
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
                'E_internal                ', &
                'L2_total                  ', &
                'L2_insideSpheriod         ', &
                'L2_outsideSpheriod        ', &
                'L2_insidePoles        '
           
           
10         format (2x,50(a25, :, 1X))
           
        else
           open (funit, file=trim(io_statsFileName), position='APPEND')
           write (funit, 11) 
11         format('# simulation restarted')
        endif
     endif
     
     write (funit, 12) simtime, gsum      ! Write the global sums to the file.
12   format (1x, 50(es25.18, :, 1x))
 
     close (funit)          ! Close the file.
    
  endif
  
  call MPI_Barrier (MPI_Comm_World, error)
  
  !=============================================================================
  
  return
end subroutine IO_writeIntegralQuantities



