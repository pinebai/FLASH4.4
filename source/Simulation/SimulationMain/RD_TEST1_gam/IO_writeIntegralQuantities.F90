!!****if* source/Simulation/SimulationMain/RD_TEST1_gam/IO_writeIntegralQuantities
!!
!!
!!  NAME
!!    IO_writeIntegralQuantities
!!
!!  SYNOPSIS
!!    call IO_writeIntegralQuantities(integer(in) :: isFirst,
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

  use IO_data, ONLY : io_restart, io_statsFileName, io_globalComm
  use Grid_interface, ONLY : Grid_getListOfBlocks, &
    Grid_getBlkIndexLimits, Grid_getBlkPtr, Grid_getSingleCellVol, &
    Grid_releaseBlkPtr, Grid_getCellCoords

   use IO_data, ONLY : io_globalMe
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

  integer, parameter ::  nGlobalSum = 17          ! Number of globally-summed quantities
  real :: gsum(nGlobalSum) !Global summed quantities
  real :: lsum(nGlobalSum) !Global summed quantities

  integer :: i, j, k
  real :: dvol             !, del(MDIM)
  real, DIMENSION(:,:,:,:), POINTER :: solnData

  integer :: point(MDIM)
  integer :: ioStat

  integer :: isize, jsize, ksize
  real,allocatable, dimension(:) :: xCenter, yCenter, zCenter

  ! Sum quantities over all locally held leaf-node blocks.
  gsum  = 0.
  lsum = 0.
  
  call Grid_getListOfBlocks(LEAF, blockList, count)
  
  do lb = 1, count
     !get the index limits of the block
     call Grid_getBlkIndexLimits(blockList(lb), blkLimits, blkLimitsGC)

     ! get a pointer to the current block of data
     call Grid_getBlkPtr(blockList(lb), solnData)

     isize = blkLimitsGC(HIGH,IAXIS)-blkLimitsGC(LOW,IAXIS)+1
     allocate(xCenter(isize))
     call Grid_getCellCoords(IAXIS,blockList(lb),&
                                 CENTER,.TRUE.,xCenter,isize)
!!     delX = xCenter(2) - xCenter(1)
#if NDIM > 1
     jsize = blkLimitsGC(HIGH,JAXIS)-blkLimitsGC(LOW,JAXIS)+1
     allocate(yCenter(jsize))
     call Grid_getCellCoords(JAXIS,blockList(lb),&
                                 CENTER,.TRUE.,yCenter,jsize)
#if NDIM == 3
     ksize = blkLimitsGC(HIGH,KAXIS)-blkLimitsGC(LOW,KAXIS)+1
     allocate(zCenter(ksize))
     call Grid_getCellCoords(KAXIS,blockList(lb),&
                                 CENTER,.TRUE.,zCenter,ksize)
#endif
#endif


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
#ifdef HE4_SPEC
#ifdef C12_SPEC
#ifdef O16_SPEC
#ifdef NE20_SPEC
#ifdef MG24_SPEC
#ifdef SI28_SPEC
#ifdef S32_SPEC
#ifdef FE54_SPEC
              lsum(2) = lsum(2) + solnData(DENS_VAR,i,j,k)*solnData(HE4_SPEC,i,j,k)*dvol
              lsum(3) = lsum(3) + solnData(DENS_VAR,i,j,k)*solnData(C12_SPEC,i,j,k)*dvol
              lsum(4) = lsum(4) + solnData(DENS_VAR,i,j,k)*solnData(O16_SPEC,i,j,k)*dvol
              lsum(5) = lsum(5) + solnData(DENS_VAR,i,j,k)*solnData(NE20_SPEC,i,j,k)*dvol
              lsum(6) = lsum(6) + solnData(DENS_VAR,i,j,k)*solnData(MG24_SPEC,i,j,k)*dvol
              lsum(7) = lsum(7) + solnData(DENS_VAR,i,j,k)*solnData(SI28_SPEC,i,j,k)*dvol
              lsum(8) = lsum(8) + solnData(DENS_VAR,i,j,k)*solnData(S32_SPEC,i,j,k)*dvol                      
	      lsum(9) = lsum(9) + solnData(DENS_VAR,i,j,k)*solnData(FE52_SPEC,i,j,k)*dvol
!              lsum(10) = lsum(10) + solnData(DENS_VAR,i,j,k)*solnData(FE54_SPEC,i,j,k)*dvol  
              lsum(10) = lsum(10) + solnData(DENS_VAR,i,j,k)*solnData(NI56_SPEC,i,j,k)*dvol
!              lsum(12) = lsum(12) + solnData(DENS_VAR,i,j,k)*solnData(CR56_SPEC,i,j,k)*dvol

            
#endif           
#endif
#endif
#endif
#endif
#endif
#endif
#endif
#endif

#ifdef DENS_VAR
#ifdef VELX_VAR      
              ! momentum
              lsum(11) = lsum(11) + solnData(DENS_VAR,i,j,k) * & 
                   &                                solnData(VELX_VAR,i,j,k)*dvol 
           
#endif
#ifdef VELY_VAR      

              lsum(12) = lsum(12) + solnData(DENS_VAR,i,j,k) * & 
                   &                                solnData(VELY_VAR,i,j,k)*dvol
           
#endif
#ifdef VELZ_VAR      
#if NDIM < 3
              lsum(13) = lsum(13) + solnData(DENS_VAR,i,j,k) * & 
                   &                                solnData(VELZ_VAR,i,j,k)*dvol*xCenter(i)
#endif
#endif

              ! total energy
#ifdef ENER_VAR
              lsum(14) = lsum(14) + solnData(ENER_VAR,i,j,k) * & 
                   &                                solnData(DENS_VAR,i,j,k)*dvol
#ifdef MAGP_VAR
              ! total plasma energy
!!$              lsum(5) = lsum(5) + (solnData(ENER_VAR,i,j,k) * & 
!!$                   &    solnData(DENS_VAR,i,j,k) + solnData(MAGP_VAR,i,j,k))*dvol

!              lsum(14) = lsum(14) + solnData(MAGP_VAR,i,j,k)*dvol
#endif
!#endif

           
#ifdef VELX_VAR      
#ifdef VELY_VAR      
#ifdef VELZ_VAR      
              ! kinetic energy
              lsum(15) = lsum(15) + 0.5*solnData(DENS_VAR,i,j,k) * & 
                   &                             (solnData(VELX_VAR,i,j,k)**2+ & 
                   &                              solnData(VELY_VAR,i,j,k)**2+ & 
                   &                              solnData(VELZ_VAR,i,j,k)**2)*dvol           

#endif
#endif
#endif


#ifdef EINT_VAR
#ifdef ENUC_VAR
              ! internal energy
              lsum(16) = lsum(16) + solnData(DENS_VAR,i,j,k) * & 
                   &                                solnData(EINT_VAR,i,j,k)*dvol
              lsum(17) = lsum(17) + solnData(DENS_VAR,i,j,k) * &
                   &                                solnData(ENUC_VAR,i,j,k)*dvol
#endif
#endif
#endif ! ifdef DENS_VAR


!#ifdef MAGP_VAR
              ! magnetic energy
!              lsum(18) = lsum(18) + solnData(MAGP_VAR,i,j,k)*dvol
#endif
!Oxygen mass
!#ifdef DENS_VAR
!#ifdef O16_SPEC
!              lsum(9) = lsum(9) + solnData(DENS_VAR,i,j,k)*dvol! * &
!                   &                      solnData(O16_SPEC,i,j,k)*dvol
            
!#endif
!#endif 
              !! Uncomment the following lines to print out all the data to the screen
!!$              write(*,*) xCenter(i), solnData(DENS_VAR,i,j,k), solnData(TEMP_VAR,i,j,k), &
!!$                   solnData(PRES_VAR,i,j,k), solnData(EINT_VAR,i,j,k), solnData(VELX_VAR,i,j,k), &
!!$                   solnData(YE_MSCALAR,i,j,k), solnData(NEUT_SPEC,i,j,k), solnData(H1_SPEC,i,j,k), &
!!$                   solnData(PROT_SPEC,i,j,k), solnData(HE3_SPEC,i,j,k), solnData(HE4_SPEC,i,j,k), &
!!$                   solnData(C12_SPEC,i,j,k), solnData(N14_SPEC,i,j,k), solnData(O16_SPEC,i,j,k), &
!!$                   solnData(NE20_SPEC,i,j,k), solnData(MG24_SPEC,i,j,k), solnData(SI28_SPEC,i,j,k), &
!!$                   solnData(S32_SPEC,i,j,k), solnData(AR36_SPEC,i,j,k), solnData(CA40_SPEC,i,j,k), &
!!$                   solnData(TI44_SPEC,i,j,k), solnData(CR48_SPEC,i,j,k), solnData(FE54_SPEC,i,j,k), &
!!$                   solnData(NI56_SPEC,i,j,k)
           enddo
        enddo
     enddo
     call Grid_releaseBlkPtr(blockList(lb), solnData)
     deallocate(xCenter)
#if NDIM > 1
     deallocate(yCenter)
#if NDIM == 3
     deallocate(zCenter)
#endif
#endif

  enddo
  
  
  ! Now the MASTER_PE sums the local contributions from all of
  ! the processors and writes the total to a file.
  
  call MPI_Reduce (lsum, gsum, nGlobalSum, FLASH_REAL, MPI_SUM, & 
       &                MASTER_PE, io_globalComm, error)
  

  if (io_globalMe  == MASTER_PE) then
     
     ! create the file from scratch if it is a not a restart simulation, 
     ! otherwise append to the end of the file
     
     !No mater what, we are opening the file. Check to see if already there
     ioStat = 0
     open(funit, file=trim(io_statsFileName), position='APPEND', status='OLD', iostat=ioStat)
     if(ioStat .NE. 0) then
        !print *, 'FILE FOUND'
        open(funit, file=trim(io_statsFileName), position='APPEND')
     endif
     
     if (isFirst .EQ. 1 .AND. (.NOT. io_restart .or. ioStat .NE. 0)) then
        
#ifndef MAGP_VAR
        write (funit, 10)               &
             '#time                     ', &
             'mass                      ', &
             'HE4_total                 ', &
             'C12_total                 ', &
             'O16_total                 ', &
             'NE20_total                ', &
             'MG24_total                ', &
             'SI28_total                ', &
             'S32_total                 ', &
	     'FE52_total		', &
             'NI56_total                ', &
             'x-momentum                ', &
             'y-momentum                ', & 
             'z-momentum                ', &
             'E_total                   ', &
             'E_kinetic                 ', &
             'E_internal                ', &
             'E_nuclear                 '
#endif
        
#ifdef MAGP_VAR
        write (funit, 10)               &
             '#time                     ', &
             'mass                      ', &
             'HE4_total                 ', &
             'C12_total                 ', &
             'O16_total                 ', &
             'NE20_total                ', &
             'MG24_total                ', &
             'SI28_total                ', &
             'S32_total                 ', &
             'FE54_total                ', &
             'x-momentum                ', &
             'y-momentum                ', & 
             'z-momentum                ', &
             'E_total                   ', &
             'E_kinetic                 ', &
             'E_internal                ', &
             'E_nuclear                 ', &
             'MagEnergy                 '
#endif
        
10         format (2x,50(a25, :, 1X))

     else if(isFirst .EQ. 1) then
        write (funit, 11) 
11      format('# simulation restarted')
     endif
     
     
     write (funit, 12) simtime, gsum      ! Write the global sums to the file.
12   format (1x, 50(es25.18, :, 1x))
 
     close (funit)          ! Close the file.
     
  endif
  
  call MPI_Barrier (io_globalComm, error)
  
  !=============================================================================
  
  return
end subroutine IO_writeIntegralQuantities



