!!****if* source/Simulation/SimulationMain/unitTest/Gravity/BHTree-cylinder/Simulation_initBlock
!!
!! NAME
!!
!!  Simulation_initBlock
!!
!! 
!! SYNOPSIS
!!
!!  Simulation_initBlock(integer :: blockId)
!!
!!
!! DESCRIPTION
!!
!!  Initializes fluid data (density, pressure, velocity, etc.) for
!!  a specified block.  This version sets up the self--gravitating 
!!  isothermal cylinder in hydrostatic equilibrium (Ostriker, J.1964,  ApJ. 140, 1056)
!!
!!
!! ARGUMENTS
!!
!!  blockId -        The number of the block to initialize
!!
!! PARAMETERS
!!
!!
!!
!!***

subroutine Simulation_initBlock (blockId)

  use Simulation_data
  use Grid_interface, ONLY : Grid_getBlkIndexLimits, Grid_getBlkPtr, Grid_releaseBlkPtr,&
    Grid_getCellCoords, Grid_putPointData
  use Multispecies_interface, ONLY : Multispecies_getSumInv
!  use PhysicalConstants_interface, ONLY : PhysicalConstants_get
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
  use Driver_interface, ONLY : Driver_abortFlash

  implicit none

#include "constants.h"
#include "Flash.h"
#include "Multispecies.h"
#include "Eos.h"
  
  integer,intent(IN) ::  blockId
  real, dimension(NSPECIES) :: massFrac

  
  integer  :: i, j, k, istat, n
  integer  :: ii, jj, kk, jhi, jlo
  real     :: vx, vy, vz, p, rho, e, ek, gpot
  real     :: xx, yy, zz, dist, distInv, frac
  real     :: gam

  real :: aconst, rho_amb, cs2, r2
  real :: phi00, phi01, rext

  real,allocatable,dimension(:) :: xCoord,yCoord,zCoord
  integer,dimension(2,MDIM) :: blkLimits,blkLimitsGC
  integer :: sizeX,sizeY,sizeZ
  integer,dimension(MDIM) :: axis
  real, dimension(:,:,:,:),pointer :: solnData
  integer,dimension(MDIM) :: startingPos
  real          :: del(MDIM)
  character(len=MAX_STRING_LENGTH) :: grav_boundary_type_x &
  & , grav_boundary_type_y, grav_boundary_type_z
  integer :: nperiodic

  logical :: gcell = .true.
     
!! ---------------------------------------------------------------------------

  call RuntimeParameters_get("grav_boundary_type_x", grav_boundary_type_x)
  call RuntimeParameters_get("grav_boundary_type_y", grav_boundary_type_y)
  call RuntimeParameters_get("grav_boundary_type_z", grav_boundary_type_z)

! check periodicity of the problem
  nperiodic=0
  if (grav_boundary_type_x.eq."periodic") nperiodic = nperiodic + 1
  if (grav_boundary_type_y.eq."periodic") nperiodic = nperiodic + 1
  if (grav_boundary_type_z.eq."periodic") nperiodic = nperiodic + 1

! In the case with unsuitable symmetry, terminate the run.
  if (nperiodic.ne.1) call Driver_abortFlash("FATAL: gravity boundaries aren't periodic in one direction")

  ! get the coordinate information for the current block from the database
  call Grid_getBlkIndexLimits(blockId,blkLimits,blkLimitsGC)
  sizeX = blkLimitsGC(HIGH,IAXIS) - blkLimitsGC(LOW,IAXIS) + 1
  allocate(xCoord(sizeX),stat=istat)
  sizeY = blkLimitsGC(HIGH,JAXIS) - blkLimitsGC(LOW,JAXIS) + 1
  allocate(yCoord(sizeY),stat=istat)
  sizeZ = blkLimitsGC(HIGH,KAXIS) - blkLimitsGC(LOW,KAXIS) + 1
  allocate(zCoord(sizeZ),stat=istat)

  if (NDIM == 3) call Grid_getCellCoords&
                      (KAXIS, blockId, CENTER, gcell, zCoord, sizeZ)
  if (NDIM >= 2) call Grid_getCellCoords&
                      (JAXIS, blockId, CENTER,gcell, yCoord, sizeY)
  call Grid_getCellCoords(IAXIS, blockId, CENTER, gcell, xCoord, sizeX)
  call Grid_getDeltas(blockId,del)

  !
  !     For each cell
  !  
#ifdef FL_NON_PERMANENT_GUARDCELLS
  call Grid_getBlkPtr(blockId,solnData)
#endif
  !print *, "coords obtained"

! generate self--gravitating isothermal cylinder
! axis of the cylinder points in the x direction

! in the case of non--isothermal EOS terminate run
  if (abs(gamma_1-1.0).gt.0.01) call Driver_abortFlash("FATAL: Attempt to build isothermal cylinder with not isothermal EOS")
! prepare useful constants
! physical quantity: _c -- inside cylinder _a -- in ambient medium outside cylinder
  cs2 = gamma_1*sim_boltz*sim_temp_c/(abar_1*sim_mH*gamma_1)
  aconst = 0.5*sim_pi*sim_newt*sim_dens_c/cs2
  rext = sqrt(sqrt(cs2*sim_dens_c/(gamma_1*sim_press_a))-1.0)/sqrt(aconst)
  phi00 = 4.0*cs2*aconst*rext**2/(1.0+aconst*rext**2)
  phi01 = 2.0*cs2*log(1.0+aconst*rext**2)-phi00*log(rext)
  sim_corrPot = abs(2.0*cs2*log(1.0+aconst*rext**2))
  rho_amb = sim_press_a*abar_2*sim_mH*gamma_2/(sim_boltz*sim_temp_a)

! fill the cells
  do k = blkLimitsGC(LOW,KAXIS), blkLimitsGC(HIGH,KAXIS)
    zz = zCoord(k)
    do j = blkLimitsGC(LOW, JAXIS), blkLimitsGC(HIGH, JAXIS)
       yy = yCoord(j)
       do i = blkLimitsGC(LOW,IAXIS), blkLimitsGC(HIGH, IAXIS)
         xx = xCoord(i)

! density profile according to eq. 57 in Ostriker 1964
         if (grav_boundary_type_x.eq."periodic") then
           r2 = yy**2+zz**2
         else if (grav_boundary_type_y.eq."periodic") then
           r2 = xx**2+zz**2
         else if (grav_boundary_type_z.eq."periodic") then
           r2 = xx**2+yy**2
         endif

         rho =  sim_dens_c/((1.0+aconst*r2)**2)
         p = rho*cs2/gamma_1
         gpot = 2.0*cs2*log(1.0+aconst*r2)

         if (p.gt.sim_press_a) then
           massFrac(1) = smallX
           massFrac(2) = 1.0-smallX
! outside cylinder is hotter homogeneous medium with lower density
         else
           rho = rho_amb
           p = sim_press_a
           gpot = phi01+phi00*log(sqrt(r2))
           massFrac(1) = 1.0-smallX
           massFrac(2) = smallX
         endif
 

! no velocity perturbation
         vx = 0.0
         vy = 0.0
         vz = 0.0

         call Multispecies_getSumInv(GAMMA, gam, massFrac(:))
         gam       = 1/gam

         ! assume gamma-law equation of state
         ek  = 0.5*(vx*vx + vy*vy + vz*vz)
         e   = p/(gam-1.)
         e   = e/rho + ek
         e   = max (e, smallP)
           
         axis(IAXIS)=i
         axis(JAXIS)=j
         axis(KAXIS)=k


#ifdef FL_NON_PERMANENT_GUARDCELLS
         solnData(SPECIES_BEGIN,i,j,k)=1.0-(NSPECIES-1)*sim_smallX
         solnData(SPECIES_BEGIN+1:SPECIES_END,i,j,k)=sim_smallX
         solnData(DENS_VAR,i,j,k)=rho
         solnData(PRES_VAR,i,j,k)=p
         solnData(ENER_VAR,i,j,k)=e
         solnData(GAME_VAR,i,j,k)=gam
         solnData(GAMC_VAR,i,j,k)=gam
         solnData(VELX_VAR,i,j,k)=vx
         solnData(VELY_VAR,i,j,k)=vy
         solnData(VELZ_VAR,i,j,k)=vz
         solnData(PANL_VAR,i,j,k)=gpot
#else

         !if there is only one species, this loop will not execute
         do n=1,NSPECIES
           call Grid_putPointData(blockId, CENTER, SPECIES_BEGIN+n-1, EXTERIOR, axis, massFrac(n))
         enddo  ! end loop over species

         call Grid_putPointData(blockId, CENTER, DENS_VAR, EXTERIOR, axis, rho)
         call Grid_putPointData(blockId, CENTER, PRES_VAR, EXTERIOR, axis, p)
         call Grid_putPointData(blockId, CENTER, ENER_VAR, EXTERIOR, axis, e)    
         call Grid_putPointData(blockId, CENTER, GAME_VAR, EXTERIOR, axis, gam)
         call Grid_putPointData(blockId, CENTER, GAMC_VAR, EXTERIOR, axis, gam)
         call Grid_putPointData(blockId, CENTER, VELX_VAR, EXTERIOR, axis, vx)
         call Grid_putPointData(blockId, CENTER, VELY_VAR, EXTERIOR, axis, vy)
         call Grid_putPointData(blockId, CENTER, VELZ_VAR, EXTERIOR, axis, vz)
         call Grid_putPointData(blockId, CENTER, PANL_VAR, EXTERIOR, axis, gpot)
#endif
        enddo
     enddo
  enddo
#ifdef FL_NON_PERMANENT_GUARDCELLS
  call Grid_releaseBlkPtr(blockID, solnData)
#endif
  deallocate(xCoord)
  deallocate(yCoord)
  deallocate(zCoord)
  return
end subroutine Simulation_initBlock
