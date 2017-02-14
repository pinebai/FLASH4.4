!!****if* source/Simulation/SimulationMain/Jeans/IO_writeIntegralQuantities
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
!!  NOTES
!!    This version also adds the analytical solution to the Jeans problem.
!!
!!
!!***

!!REORDER(4):solnData

subroutine IO_writeIntegralQuantities ( isFirst, simTime)

  use IO_data, ONLY : io_restart, io_statsFileName
  use Grid_interface, ONLY : Grid_getListOfBlocks, &
    Grid_getBlkIndexLimits, Grid_getBlkPtr, Grid_getSingleCellVol, &
    Grid_releaseBlkPtr

  use Simulation_data, ONLY : sim_rho0, sim_p0, sim_lambdaX, sim_lambdaY, sim_lambdaZ, &
       sim_A, sim_gamma, sim_U0
!!  sim_p0       Initial ambient pressure
!!  sim_rho0     Initial ambient density
!!  sim_lambdaX  Perturbation X-wavelength
!!  sim_lambdaY  Perturbation Y-wavelength
!!  sim_lambdaZ  Perturbation Z-wavelength
!!  sim_A        Perturbation amplitude


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

  integer, parameter ::  nGlobalSum = 11          ! Number of globally-summed quantities
  real :: gsum(nGlobalSum) !Global summed quantities
  real :: lsum(nGlobalSum) !Global summed quantities

  integer :: i, j, k
  real :: dvol             !, del(MDIM)
  real, DIMENSION(:,:,:,:), POINTER :: solnData
  real :: trueKinetic, trueInternal, truePotential, G, L, c0, omega, kWave, kJeans

  integer :: point(MDIM)

#ifdef MAGP_VAR
#define MAGP_INDEX 8
#ifdef GPOT_VAR
#define GPOT_INDEX 9
#endif
#endif

#ifndef MAGP_VAR
#ifdef GPOT_VAR
#define GPOT_INDEX 8
#endif
#endif


  ! Constants for the analytical solution 
  G = 6.67428E-8  !! Newton's constant
  L = 0.5 * sqrt ((PI * sim_gamma * sim_rho0) / (G * (sim_p0**2)))
  c0 = sqrt ((sim_gamma * sim_rho0) / sim_p0) ! unperturbed adiabatic sound speed
  kJeans = sqrt(4*PI*G*sim_p0) / c0
  kWave = 2*PI/sim_lambdaX   !! 10.984  Only a one-dimensional oscillation in X
  omega = sqrt (((c0**2) * (kWave**2)) - (4.0*PI*G*sim_p0)) 


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
              ! total plasma energy
#ifdef MAGP_VAR
              lsum(5) = lsum(5) + (solnData(MAGP_VAR,i,j,k))*dvol
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

#ifdef MAGP_VAR
              ! magnetic energy
              lsum(MAGP_INDEX) = lsum(MAGP_INDEX) + solnData(MAGP_VAR,i,j,k)*dvol
#endif

#ifdef GPOT_VAR
              ! gravitational potential energy
              lsum(GPOT_INDEX) = lsum(GPOT_INDEX) + 0.5 * &
                   &                                (solnData(DENS_VAR,i,j,k) &
                   &                                - sim_rho0)* & 
                   &                                solnData(GPOT_VAR,i,j,k) * &
                   &                                dvol
              lsum(5) = lsum(5) + 0.5*solnData(DENS_VAR,i,j,k) * & 
                   &                                solnData(GPOT_VAR,i,j,k) * &
                   &                                dvol
#endif

!! Analytical solutions, taken from the file JeansAnalytical.F90
            
              !! Kinetic energy T (versus calculated at lsum(6))
              !!   T = (p0 * (delta**2) * (omega**2) * (L**2) * (1.0 - cos(2.0*omega*time))) / &
              !!          (8.0 * (k**2))
              lsum(9) = (sim_p0*(sim_A**2) * (omega**2) * (L**2)) / (8.0*(kWave**2)) * &
                   (1.0 - cos(2.0*omega*simTime)) 


              !! Internal energy U(t) - U0
              !! U = -0.125 * p0 * (c0**2) * (delta**2) * (L**2) * (1.0 - cos(2.0*omega*time))
              lsum(10) = -0.125 * sim_p0 * (c0**2) * (sim_A**2) * (L**2) * (1.0 - cos(2.0*omega*simTime))
              !  Subtract off initial internal energy
              if (isFirst .eq. 1) then
                 sim_U0 = lsum(10)
                !DEV: commented out to prevent flooding test-suite with 65-meg ttext files --PR 
                 !print *, 'First time through, U0 = ',sim_U0
              endif
              !DEV: See above --PR
              !print *, 'Second time through, U0 = ',sim_U0
              lsum(10) = lsum(10) - sim_U0

              !! Potential energy W
              !! W = -((pi * G * (p0**2) * (delta**2) * (L**2)) / (2.0 * (k**2))) * (1.0 + cos(2.0*omega*time))
              lsum(11) = - ((PI * G * (sim_p0**2) * (sim_A**2) * (L**2)) / &
                   (2.0 * (kWave**2))) * (1.0 + cos(2.0*omega*simTime))
           enddo
        enddo
     enddo
     call Grid_releaseBlkPtr(blockList(lb), solnData)

  enddo
  

  
  ! Now the MASTER_PE sums the local contributions from all of
  ! the processors and writes the total to a file.
  
  call MPI_Reduce (lsum, gsum, nGlobalSum, FLASH_REAL, MPI_SUM, & 
       &                MASTER_PE, MPI_COMM_WORLD, error)
  

  if (io_globalMe  == MASTER_PE) then
     
     ! create the file from scratch if it is a not a restart simulation, 
     ! otherwise append to the end of the file
     if (isfirst == 0) then
        open (funit, file=trim(io_statsFileName), position='APPEND')
     else 
        if (.NOT. io_restart) then
           open (funit, file=trim(io_statsFileName)) 
#ifndef MAGP_VAR
#ifdef GPOT_VAR
           write (funit, 10)               &
                '#time                     ', &
                'mass                      ', &
                'x-momentum                ', &
                'y-momentum                ', & 
                'z-momentum                ', &
                'E_total                   ', &
                'E_kinetic                 ', &
                'E_internal                ', &
                'E_gravity                 ', &
                'E_kinetic analytical      ', &
                'E_internal analytical     ', &
                'E_gravity analytical      '
#else
           write (funit, 10)               &
                '#time                     ', &
                'mass                      ', &
                'x-momentum                ', &
                'y-momentum                ', & 
                'z-momentum                ', &
                'E_total                   ', &
                'E_kinetic                 ', &
                'E_internal                ', &
                'E_kinetic analytical      ', &
                'E_internal analytical     ', &
                'E_gravity analytical      '
#endif

#ifdef MAGP_VAR

           write (funit, 10)               &
                '#time                     ', &
                'mass                      ', &
                'x-momentum                ', &
                'y-momentum                ', & 
                'z-momentum                ', &
                'E_total                   ', &
                'E_kinetic                 ', &
                'E_internal                ', &
                'E_magnetic                ', &
                'E_gravity                 ', &
                'E_kinetic analytical      ', &
                'E_internal analytical     ', &
                'E_gravity analytical      '
#else
           write (funit, 10)               &
                '#time                     ', &
                'mass                      ', &
                'x-momentum                ', &
                'y-momentum                ', & 
                'z-momentum                ', &
                'E_total                   ', &
                'E_kinetic                 ', &
                'E_internal                ', &
                'E_magnetic                ', &
                'E_kinetic analytical      ', &
                'E_internal analytical     ', &
                'E_gravity analytical      '
#endif
#endif

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



