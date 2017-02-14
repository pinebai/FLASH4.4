!!****if* source/Simulation/SimulationMain/unitTest/SinkMomTest/IO_writeIntegralQuantities
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
    Grid_releaseBlkPtr
  use Particles_interface, ONLY : Particles_sinkSumAttributes

   use IO_data, ONLY : io_globalMe, io_writeMscalarIntegrals
   use Simulation_data, ONLY : sim_testInitialized, &
        sim_testRefVals, &
        sim_testLastVals, &
        sim_testWorstVals, &
        sim_testWorstWhen
  implicit none

#include "Flash_mpi.h"
#include "constants.h"
#include "Flash.h"
  
  
  real, intent(in) :: simTime

  integer, intent(in) :: isFirst

  integer :: lb, count
  
  integer :: funit = 99
  integer :: error
  integer :: nGlobalSumUsed, iSum
  
  character (len=MAX_STRING_LENGTH), save :: fname 
  
  integer :: blockList(MAXBLOCKS)

  integer :: blkLimits(HIGH, MDIM), blkLimitsGC(HIGH, MDIM)

#ifdef MAGP_VAR
  integer, parameter ::  nGlobalSumProp = 8              ! Number of globally-summed regular quantities
#else
  integer, parameter ::  nGlobalSumProp = 7              ! Number of globally-summed regular quantities
#endif
  integer, parameter ::  nGlobalSum = nGlobalSumProp + NMASS_SCALARS ! Number of globally-summed quantities
  real :: gsum(nGlobalSum) !Global summed quantities
  real :: lsum(nGlobalSum) !Global summed quantities


  real,   dimension(1) :: sinkPartSums0
  real,   dimension(3) :: sinkPartSums
  integer,dimension(1),save :: sinkPartAttribs0=(/MASS_PART_PROP/)
  integer,dimension(3),save :: sinkPartAttribs=(/VELX_PART_PROP,VELY_PART_PROP,VELZ_PART_PROP/)
  real :: gsumExtra(1+3) !Extra summed quantities, includes sink particle contributions

  integer :: ivar
  integer :: i, j, k
  real :: dvol             !, del(MDIM)
  real, DIMENSION(:,:,:,:), POINTER :: solnData

  integer :: point(MDIM)
  integer :: ioStat

  if (io_writeMscalarIntegrals) then
     nGlobalSumUsed = nGlobalSum
  else
     nGlobalSumUsed = nGlobalSumProp
  end if

  ! Sum quantities over all locally held leaf-node blocks.
  gsum(1:nGlobalSumUsed) = 0.
  lsum(1:nGlobalSumUsed) = 0.
  
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
#ifdef MAGP_VAR
              ! total plasma energy
!!$              lsum(5) = lsum(5) + (solnData(ENER_VAR,i,j,k) * & 
!!$                   &    solnData(DENS_VAR,i,j,k) + solnData(MAGP_VAR,i,j,k))*dvol

              lsum(5) = lsum(5) + solnData(MAGP_VAR,i,j,k)*dvol
#endif
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
#endif

#ifdef MAGP_VAR
              ! magnetic energy
              lsum(8) = lsum(8) + solnData(MAGP_VAR,i,j,k)*dvol
#endif

#ifdef DENS_VAR
              if (io_writeMscalarIntegrals) then
                 iSum = nGlobalSumProp
!!$                 do ivar=MASS_SCALARS_BEGIN,MASS_SCALARS_END
                    lsum(iSum+1:iSum+NMASS_SCALARS) = &
                         lsum(iSum+1:iSum+NMASS_SCALARS) + &
                           solnData(DENS_VAR,i,j,k) * & 
                           solnData(MASS_SCALARS_BEGIN: &
                                    MASS_SCALARS_END,i,j,k)*dvol
!!$                 end do
              end if
#endif
           enddo
        enddo
     enddo
     call Grid_releaseBlkPtr(blockList(lb), solnData)

  enddo
  

  
  ! Now the MASTER_PE sums the local contributions from all of
  ! the processors and writes the total to a file.
  
  call MPI_Reduce (lsum, gsum, nGlobalSumUsed, FLASH_REAL, MPI_SUM, & 
       &                MASTER_PE, io_globalComm, error)
  

  call Particles_sinkSumAttributes(sinkPartSums0,sinkPartAttribs0)
  call Particles_sinkSumAttributes(sinkPartSums,sinkPartAttribs,MASS_PART_PROP)

  if (io_globalMe  == MASTER_PE) then
     

     gsumExtra(1) = gsum(1) + sinkPartSums0(1)
     gsumExtra(2) = gsum(2) + sinkPartSums(1)
     gsumExtra(3) = gsum(3) + sinkPartSums(2)
     gsumExtra(4) = gsum(4) + sinkPartSums(3)

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
             'mass TOTAL                ', &
             'x-momentum TOTAL          ', &
             'y-momentum TOTAL          ', & 
             'z-momentum TOTAL          ', &
             'mass (fluid)              ', &
             'x-momentum (fluid)        ', &
             'y-momentum (fluid)        ', & 
             'z-momentum (fluid ...)    ', &
             'E_total                   ', &
             'E_kinetic                 ', &
             'E_internal                ', &
             (msName(ivar),ivar=MASS_SCALARS_BEGIN,&
              min(MASS_SCALARS_END,&
                  MASS_SCALARS_BEGIN+nGlobalSumUsed-nGlobalSumProp-1))

#else
        
        write (funit, 10)               &
             '#time                     ', &
             'mass TOTAL                ', &
             'x-momentum TOTAL          ', &
             'y-momentum TOTAL          ', & 
             'z-momentum TOTAL          ', &
             'mass (fluid)              ', &
             'x-momentum (fluid)        ', &
             'y-momentum (fluid)        ', & 
             'z-momentum (fluid ...)    ', &
             'E_total                   ', &
             'E_kinetic                 ', &
             'E_internal                ', &
             'MagEnergy                 ', &
             (msName(ivar),ivar=MASS_SCALARS_BEGIN,&
              min(MASS_SCALARS_END,&
                  MASS_SCALARS_BEGIN+nGlobalSumUsed-nGlobalSumProp-1))
#endif
        
10         format (2x,50(a25, :, 1X))

     else if(isFirst .EQ. 1) then
        write (funit, 11) 
11      format('# simulation restarted')
     endif
     
     ! Write the global sums to the file.
     write (funit, 12) simtime, gsumExtra(1:4), gsum(1:nGlobalSumUsed)

12   format (1x, 50(es25.18, :, 1x))
 
     close (funit)          ! Close the file.
     
     if (.NOT.sim_testInitialized) then
        sim_testRefVals(1:4) = gsumExtra(1:4)
        sim_testWorstVals(1:4) = gsumExtra(1:4)
        sim_testWorstWhen(1:4) = simTime
        sim_testInitialized = .TRUE.
     end if
     sim_testLastVals(1:4) = gsumExtra(1:4)
     do i=1,4
        if (abs(sim_testLastVals(i)-sim_testRefVals(i)) >     &
            abs(sim_testWorstVals(i)-sim_testRefVals(i))) then
            sim_testWorstVals(i) = sim_testLastVals(i)
            sim_testWorstWhen(i) = simTime
         end if
     end do


  endif


  
  call MPI_Barrier (io_globalComm, error)
  
  !=============================================================================
  
  return

  contains
    character(len=25) function msName(ivar)
      integer,intent(in) :: ivar
      character(len=25) :: str
      call Simulation_mapIntToStr(ivar,str,MAPBLOCK_UNK)
      msName = str
    end function msName
end subroutine IO_writeIntegralQuantities



