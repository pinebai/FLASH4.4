!!****if* source/Simulation/SimulationMain/Plasma/IO_writeIntegralQuantities
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

subroutine IO_writeIntegralQuantities ( isFirst, simTime)
  use Simulation_data, ONLY : sim_pmass_1, sim_pmass_2, &
       sim_esw, sim_bsw, sim_usw
  use Particles_data, only: particles, pt_numLocal
  use IO_data, only : io_restart, io_statsFileName
  use Grid_interface, only : Grid_getListOfBlocks, &
    Grid_getBlkIndexLimits, Grid_getBlkPtr, Grid_getSingleCellVol, &
    Grid_releaseBlkPtr, Grid_getBlkPhysicalSize, Grid_getBlkIndexLimits
   use IO_data, ONLY : io_globalMe
  implicit none

#include "Flash_mpi.h"
#include "constants.h"
#include "Flash.h"
#include "Particles.h"

  
  real, intent(in) :: simTime

  integer, intent(in) :: isFirst

  ! vacuum permittivity [F/m] (electric constant) 
  real, parameter :: e0 = 8.854187817e-12
  ! permeability of vacuum [Vs/Am] (magnetic constant) 
  real, parameter :: mu0 = 4.0e-7*3.141592653589793238462643383279503
  real, parameter :: amu = 1.66054e-27   ! Atomic mass unit [kg]
  
  integer :: lb, count
  
  integer :: funit = 99
  integer :: error
  
  character (len=MAX_STRING_LENGTH), save :: fname 
  
  integer :: blockList(MAXBLOCKS), blockCount, localSize(MDIM)
  integer :: blkLimits(HIGH, MDIM), blkLimitsGC(HIGH, MDIM)
  real :: blockSize(MDIM), blockCenter(MDIM)

  logical, save :: first_call = .true.

  real :: h(MDIM), dv, b(MDIM), e(MDIM), divb
  real :: vb(3,2), n(2), nprot 
  real :: v(MDIM), nbsw(MDIM), vpara, vperp(MDIM)

  integer, parameter ::  nGlobalSum = 15 ! No. of globally-summed quantities
  real :: gsum(nGlobalSum) !Global summed quantities
  real :: lsum(nGlobalSum) !Global summed quantities

  real, DIMENSION(:,:,:,:), POINTER :: u
  integer :: i, j, k, s

  integer :: point(MDIM)
  integer :: ioStat


  ! Sum quantities over all locally held leaf-node blocks.
  gsum  = 0.
  lsum = 0.
  
  call Grid_getListOfBlocks(LEAF, blockList, count)
  
  do lb = 1, count
     !get the index limits of the block
     call Grid_getBlkIndexLimits(blockList(lb), blkLimits, blkLimitsGC)

     ! get a pointer to the current block of data
     call Grid_getBlkPtr(blockList(lb), u)
     call Grid_getBlkPhysicalSize(blockList(lb), blockSize)
     call Grid_getBlkIndexLimits(blockList(lb), &
          blkLimits, blkLimitsGC, CENTER)
     localSize=blkLimits(HIGH,:)-blkLimits(LOW,:)+1   ! NXB, NYB and NZB

     h = blockSize/localSize  ! cell size
     dv = product(h(1:NDIM))  ! cell volume [m^3]

     ! Sum contributions from the indicated blkLimits of cells.
     do k = blkLimits(LOW,KAXIS), blkLimits(HIGH,KAXIS)
        do j = blkLimits(LOW,JAXIS), blkLimits(HIGH,JAXIS)
           do i = blkLimits(LOW,IAXIS), blkLimits(HIGH,IAXIS)
              
              point(IAXIS) = i
              point(JAXIS) = j
              point(KAXIS) = k

              ! Get the cell volume for a single cell
              !call Grid_getSingleCellVol(blockList(lb), EXTERIOR, point, dv)
              
              ! Accumulate energy of the magnetic field. 
              ! -sim_bsw() for energy of the fluctuating field
              b(1) = u(GRBX_VAR,i,j,k) !-sim_bsw(1)
              b(2) = u(GRBY_VAR,i,j,k) !-sim_bsw(2)
              b(3) = u(GRBZ_VAR,i,j,k) !-sim_bsw(3)
              lsum(1) = lsum(1) + dv*(sum(b*b))*0.5/mu0 
              ! Accumulate energy of the electric field. 
              ! -sim_esw() for energy of the fluctuating field
              e(1) = u(GREX_VAR,i,j,k) !-sim_esw(1) 
              e(2) = u(GREY_VAR,i,j,k) !-sim_esw(2)
              e(3) = u(GREZ_VAR,i,j,k) !-sim_esw(3)
              lsum(2) = lsum(2) + dv*(sum(e*e))*0.5*e0 
              ! Accumulate div(b)
              ! for ndim==1 div(b)=0 by definition
              divb = (u(GRBX_VAR,i+1,j,k)-u(GRBX_VAR,i-1,j,k))/(2*h(1))
              if (NDIM > 1) then
                 divb = divb &
                      + (u(GRBY_VAR,i,j+1,k)-u(GRBY_VAR,i,j-1,k))/(2*h(2))
              end if
              if (NDIM == 3) then
                 ! should use same stencil as in rot for div
                 divb = divb &
                      + (u(GRBZ_VAR,i,j,k+1)-u(GRBZ_VAR,i,j,k-1))/(2*h(3))
              end if
              lsum(15) = lsum(15) + abs(divb)
           enddo
        enddo
     enddo
     call Grid_releaseBlkPtr(blockList(lb), u)
  enddo

  ! Now compute particle quantities
  ! compute bulk velocity for species 1 and 2 for thermal energy below
  ! note that this is only the average on this processor...
  vb = 0.0
  n = 0.0
  do i = 1, pt_numLocal
     s = particles(SPECIE_PART_PROP, i)
     ! number of real protons
     if (s == 1) then 
        nprot = particles(MASS_PART_PROP, i)/(amu*sim_pmass_1) 
     else if (s == 2) then 
        nprot = particles(MASS_PART_PROP, i)/(amu*sim_pmass_2) 
     end if     
     if (s < 3) then 
        vb(1,s) = vb(1,s) + nprot*particles(VELX_PART_PROP, i)
        vb(2,s) = vb(2,s) + nprot*particles(VELY_PART_PROP, i) 
        vb(3,s) = vb(3,s) + nprot*particles(VELZ_PART_PROP, i) 
        n(s) = n(s) + nprot
     end if
  end do
  vb(:,1) = vb(:,1)/n(1)
  vb(:,2) = vb(:,2)/n(2)

  do i = 1, pt_numLocal
     ! s==1 for core, and s==2 for beam protons
     s = particles(SPECIE_PART_PROP, i)
     v = (/ particles(VELX_PART_PROP, i), particles(VELY_PART_PROP, i), &
          particles(VELZ_PART_PROP, i) /)
     ! number of real protons
     if (s == 1) then 
        nprot = particles(MASS_PART_PROP, i)/(amu*sim_pmass_1) 
     else if (s == 2) then 
        nprot = particles(MASS_PART_PROP, i)/(amu*sim_pmass_2) 
     end if
     if (s < 3) then 
        ! Accumulate total kinetic energy 
        lsum(2+s) = lsum(2+s) + 0.5*particles(MASS_PART_PROP, i)*sum(v*v)
        ! Accumulate number of meta particles
        lsum(4+s) = lsum(4+s) + 1.0
        ! Accumulate real number of protons
        lsum(6+s) = lsum(6+s) + nprot

        nbsw = sim_bsw/sqrt(sum(sim_bsw*sim_bsw))   ! normalized magnetic field vector
        vpara = dot_product(v, nbsw)    ! v.b, field parallel velocity
        vperp = v-vpara*nbsw            ! perpendicular velocity
        ! Accumulate parallel velocity
        lsum(8+s) = lsum(8+s) + nprot*vpara

        v = v - vb(:,s)              ! thermal velocity
        vpara = dot_product(v, nbsw) ! v.b, field parallel thermal velocity
        vperp = v-vpara*nbsw         ! perpendicular thermal velocity
        ! Accumulate total parallel thermal energy 
        lsum(10+s) = lsum(10+s) + &
             0.5*particles(MASS_PART_PROP, i)*vpara*vpara
        ! Accumulate total perpendicular thermal energy 
        lsum(12+s) = lsum(12+s) + &
             0.5*particles(MASS_PART_PROP, i)*sum(vperp*vperp)
     end if
  end do
  

  
  ! Now the MASTER_PE sums the local contributions from all of
  ! the processors and writes the total to a file.
  
  call MPI_Reduce (lsum, gsum, nGlobalSum, FLASH_REAL, MPI_SUM, & 
       &                MASTER_PE, MPI_COMM_WORLD, error)
  

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
        write (funit, 10)               &
             '#time                     ', &
             'magnetic_energy           ', &
             'electric_energy           ', &
             'kinetic_energy_1          ', & 
             'kinetic_energy_2          ', &
             'meta_1                    ', &
             'meta_2                    ', &
             'protons_1                 ', &
             'protons_2                 ', &
             'xvelsum_1                 ', &
             'xvelsum_2                 ', &
             'thermal_para_1            ', &
             'thermal_para_2            ', &
             'thermal_perp_1            ', &
             'thermal_perp_2            ', &
             'divbsum                   '
10         format (2x,50(a25, :, 1X))

     else if(isFirst .EQ. 1) then
        write (funit, 11) 
11      format('# simulation restarted')
     endif
     
     
     write (funit, 12) simtime, gsum      ! Write the global sums to the file.
12   format (1x, 50(es25.18, :, 1x))
 
     close (funit)          ! Close the file.
     
  endif
  
  call MPI_Barrier (MPI_Comm_World, error)
  
  !=============================================================================
  
  return
end subroutine IO_writeIntegralQuantities
