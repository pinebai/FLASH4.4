!!****if* source/Simulation/SimulationMain/unitTest/Gravity/BHTree-layer/Simulation_initBlock
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
!!  a specified block.  This version sets up the spherical expanding shell.
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
  
  implicit none

#include "constants.h"
#include "Flash.h"
#include "Multispecies.h"
#include "Eos.h"
  
  integer,intent(IN) ::  blockId

  
  integer  :: i, j, k, istat, n
  integer  :: ii, jj, kk, jhi, jlo
  real     :: vx, vy, vz, p, rho, e, ek, gpot
  real     :: xx, yy, zz, dist, distInv, frac

  real,allocatable,dimension(:) :: xCoord,yCoord,zCoord
  integer,dimension(2,MDIM) :: blkLimits,blkLimitsGC
  integer :: sizeX,sizeY,sizeZ
  integer,dimension(MDIM) :: axis
  real, dimension(:,:,:,:),pointer :: solnData
  integer,dimension(MDIM) :: startingPos
  real          :: del(MDIM)

  logical :: gcell = .true.

     
!! ---------------------------------------------------------------------------

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


  do k = blkLimitsGC(LOW,KAXIS), blkLimitsGC(HIGH,KAXIS)
    zz = zCoord(k)
    do j = blkLimitsGC(LOW, JAXIS), blkLimitsGC(HIGH, JAXIS)
       yy = yCoord(j)
       do i = blkLimitsGC(LOW,IAXIS), blkLimitsGC(HIGH, IAXIS)
         xx = xCoord(i)
             
         if (sim_dir == 1) then
           dist = xx
         else 
           if (sim_dir == 2) then
             dist = yy
           else
             dist = zz
           endif
         endif
         call sim_find (sim_zProf, sim_nProfile, abs(dist), jlo)


         if (jlo .eq. 0) then
            jlo = 1
            jhi = 1
            frac = 0.
         else if (jlo .eq. sim_nProfile) then
            jlo = sim_nProfile
            jhi = sim_nProfile
            frac = 0.
         else
            jhi = jlo + 1
            frac = (abs(dist) - sim_zProf(jlo)) / & 
                 (sim_zProf(jhi)-sim_zProf(jlo))
         endif
                 
         p    = max(sim_pProf(jlo) + frac*(sim_pProf(jhi)  - sim_pProf(jlo)), smallp)
         rho  = max(sim_rhoProf(jlo) + frac*(sim_rhoProf(jhi)- sim_rhoProf(jlo)), smlrho)
         gpot = sim_PhiProf(jlo) + frac*(sim_PhiProf(jhi)- sim_PhiProf(jlo))
         if (sim_dir == 1) then
           vx   = sim_vProf(jlo) + frac*(sim_vProf(jhi)  - sim_vProf(jlo))
           vy   = 0.0
           vz   = 0.0
         else 
           if (sim_dir == 2) then
             vx   = 0.0
             vy   = sim_vProf(jlo) + frac*(sim_vProf(jhi)  - sim_vProf(jlo))
             vz   = 0.0
           else
             vx   = 0.0
             vy   = 0.0
             vz   = sim_vProf(jlo) + frac*(sim_vProf(jhi)  - sim_vProf(jlo))
           endif
         endif

         ! assume gamma-law equation of state
         ek  = 0.5*(vx*vx + vy*vy + vz*vz)
         e   = p/(sim_gamma_1-1.)
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
         solnData(GAME_VAR,i,j,k)=sim_gamma_1
         solnData(GAMC_VAR,i,j,k)=sim_gamma_1
         solnData(VELX_VAR,i,j,k)=vx
         solnData(VELY_VAR,i,j,k)=vy
         solnData(VELZ_VAR,i,j,k)=vz
         solnData(PANL_VAR,i,j,k)=gpot
#else


         call Grid_putPointData(blockId, CENTER, DENS_VAR, EXTERIOR, axis, rho)
         call Grid_putPointData(blockId, CENTER, PRES_VAR, EXTERIOR, axis, p)
         call Grid_putPointData(blockId, CENTER, ENER_VAR, EXTERIOR, axis, e)    
         call Grid_putPointData(blockId, CENTER, GAME_VAR, EXTERIOR, axis, sim_gamma_1)
         call Grid_putPointData(blockId, CENTER, GAMC_VAR, EXTERIOR, axis, sim_gamma_1)
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






!******************************************************************************

!  Routine:     sim_find()

!  Description: Given a monotonically increasing table x(N) and a test value
!               x0, return the index i of the largest table value less than
!               or equal to x0 (or 0 if x0 < x(1)).  Use binary search.

subroutine sim_find (x, Narr, x0, i)

  implicit none

! Arguments, LBR guessed intent on these
  integer, intent(IN) :: Narr
  integer, intent(OUT):: i
  real, intent(IN)    :: x(Narr), x0

! local variables
  integer  il, ir, im

  if (x0 .lt. x(1)) then

     i = 0

  elseif (x0 .gt. x(Narr)) then

     i = Narr

  else

     il = 1
     ir = Narr
10   if (ir .eq. il+1) goto 20
     im = (il + ir) / 2
     if (x(im) .gt. x0) then
        ir = im
     else
        il = im
     endif
     goto 10
20   i = il

  endif

  return
end subroutine sim_find




