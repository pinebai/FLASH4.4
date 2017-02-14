!!****if* source/Simulation/SimulationMain/unitTest/Gravity/BHTree/Simulation_initBlock
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
  real, dimension(NSPECIES) :: massFrac

  
  integer  :: i, j, k, istat, n
  integer  :: ii, jj, kk, jhi, jlo
  real     :: vx, vy, vz, p, rho, e, ek, gpot
  real     :: xx, yy, zz, gam, dist, distInv, frac
  real :: xdist, ydist, zdist
  real :: theta, phi, cph, cth, sth, sph, sthInv
  real :: Ylm, dYdt, dYdp, rr, vel
  real :: insubzones, insubzm1, inszd
  real ::  sum_rho, sum_p, sum_vx, sum_vy, sum_vz, sum_x, sum_gpot

  real,allocatable,dimension(:) :: xCoord,yCoord,zCoord
  integer,dimension(2,MDIM) :: blkLimits,blkLimitsGC
  integer :: sizeX,sizeY,sizeZ
  integer,dimension(MDIM) :: axis
  real, dimension(:,:,:,:),pointer :: solnData
  integer,dimension(MDIM) :: startingPos
  real          :: del(MDIM)


  logical :: gcell = .true.

     
!! ---------------------------------------------------------------------------

  if (sim_nSubZones .le. 1) sim_nSubZones = 2
  insubzones = 1./real(sim_nSubZones)
  insubzm1   = 1./real(sim_nSubZones-1)
  inszd = insubzones**MDIM
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

#if NSPECIES > 0
     massFrac(1) = 1.0
     massFrac(2:NSPECIES) = 0.0
#endif

  do k = blkLimitsGC(LOW,KAXIS), blkLimitsGC(HIGH,KAXIS)
    zz = zCoord(k)
    do j = blkLimitsGC(LOW, JAXIS), blkLimitsGC(HIGH, JAXIS)
       yy = yCoord(j)
       do i = blkLimitsGC(LOW,IAXIS), blkLimitsGC(HIGH, IAXIS)
         xx = xCoord(i)
             
         sum_rho  = 0.
         sum_p    = 0.
         sum_vx   = 0.
         sum_vy   = 0.
         sum_vz   = 0.
         sum_x    = 0.
         sum_gpot = 0.

!    
!        Break the cell into nsubzones^ndim sub-zones, and look up the
!        appropriate quantities along the 1d profile for that subzone.  
!    
!        Have the final values for the zone be equal to the average of
!        the subzone values.
!  
     
         do kk = 0, (sim_nSubZones-1)
           zz    = zCoord(k) + (kk*insubzm1-.5)*del(KAXIS)
           zdist = (zz - sim_zCenter)
     
           do jj = 0, (sim_nSubZones-1)
             yy    = yCoord(j) + (jj*insubzm1-.5)*del(JAXIS)
             ydist = (yy - sim_yCenter)
     
             do ii = 0, (sim_nSubZones-1)
               xx    = xCoord(i) + (ii*insubzm1-.5)*del(IAXIS)
               xdist = xx - sim_xCenter
             
               if (NDIM .eq. 1) then
                 dist = abs(xdist)
               else
                 dist    = sqrt( xdist**2 + ydist**2 + zdist**2 )
               endif
               distinv = 1. / max( dist, 1.d-10 )
               call sim_find (sim_rProf, sim_nProfile, dist, jlo)

!
!  a point  `dist' is frac-way between jlo and jhi.   We do a
!  linear ierpolation of the quantities at jlo and jhi and sum those.
! 


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
                  frac = (dist - sim_rProf(jlo)) / & 
                       (sim_rProf(jhi)-sim_rProf(jlo))
               endif
                 
! 
!   Now total these quantities.   Note that  v is a radial velocity; 
!   we multiply by the tangents of the appropriate angles to get
!   the projections in the x, y and z directions.
!
               sum_p = sum_p +  & 
                     & sim_pProf(jlo) + frac*(sim_pProf(jhi)  - sim_pProf(jlo))
       
               sum_rho = sum_rho + & 
                       & sim_rhoProf(jlo) + frac*(sim_rhoProf(jhi)- sim_rhoProf(jlo))
       
               sum_x   = sum_x   + & 
                       & sim_xProf(jlo) + frac*(sim_xProf(jhi)- sim_xProf(jlo))

               sum_gpot = sum_gpot + & 
                       & sim_PhiProf(jlo) + frac*(sim_PhiProf(jhi)- sim_PhiProf(jlo))
       
               vel = sim_vProf(jlo) + frac*(sim_vProf(jhi)  - sim_vProf(jlo))
       
               sum_vx  = sum_vx  + vel*xdist*distinv
               sum_vy  = sum_vy  + vel*ydist*distinv
               sum_vz  = sum_vz  + vel*zdist*distinv
       
             enddo
           enddo
         enddo

         rho = max(sum_rho * inszd, smlrho)
         p   = max(sum_p   * inszd, smallp)
         gpot = sum_gpot * inszd
         vx  = sim_vx + sum_vx  * inszd
         vy  = sim_vy + sum_vy  * inszd
         vz  = sim_vz + sum_vz  * inszd
         massFrac(1) = max(smallX, sum_x * inszd)
         massFrac(2) = 1.0 - massFrac(1)

         ! include perturbations of density
         if (massFrac(1) .gt. 0.0) then 
           if (sim_pertType .eq. 1) then
             call RANDOM_NUMBER(rr)
             rho = rho * (1.0 + sim_pertamp*(2.0*rr-1.0))
             call RANDOM_NUMBER(rr)
             vx = vx + sim_velamp * (2.0*rr-1.0)
             call RANDOM_NUMBER(rr)
             vy = vy + sim_velamp * (2.0*rr-1.0)
             call RANDOM_NUMBER(rr)
             vz = vz + sim_velamp * (2.0*rr-1.0)
           endif
           if (sim_pertType .eq. 2) then
             xDist = (xCoord(i) - sim_xCenter)
             yDist = (yCoord(j) - sim_yCenter)
             zDist = (zCoord(k) - sim_zCenter)
             dist = sqrt(xDist**2 + yDist**2 + zDist**2)
             distInv = 1.0 / (dist + tiny(0.0))
             theta = acos(zDist/dist)
             phi = atan2(yDist, xDist)
             sth = sin(theta)
             cth = cos(theta)
             sph = sin(phi)
             cph = cos(phi)
             sthInv = 1.0 / (sth + tiny(0.0))
             call sim_spharm(sim_spharm_l1, sim_spharm_m1, theta, phi, Ylm, dYdt, dYdp)
             rho = rho * (1.0 + sim_pertamp / plgndr_norm * Ylm)
             !print *, "vamp, norm, dYdt, dYdp = ", sim_velamp, plgndr_norm, dYdt, dYdp
             ! normalization of grad(Ylm): use norm const of Ylm (plgndr_norm) and divide by l
             vx = vx + sim_velamp * (dYdt*cth*cph - sthInv*dYdp*sph) / plgndr_norm / sim_spharm_l1
             vy = vy + sim_velamp * (dYdt*cth*sph + sthInv*dYdp*cph) / plgndr_norm / sim_spharm_l1
             vz = vz + sim_velamp * (-dYdt*sth    + 0.0            ) / plgndr_norm / sim_spharm_l1
           endif
         endif
         
!
!  assume gamma-law equation of state
!

      
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
#if NSPECIES > 0
         solnData(SPECIES_BEGIN,i,j,k)=1.0-(NSPECIES-1)*sim_smallX
         solnData(SPECIES_BEGIN+1:SPECIES_END,i,j,k)=sim_smallX
#endif
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




