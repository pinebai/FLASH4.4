!!****if* source/physics/Hydro/HydroMain/unsplit/MHD_StaggeredMesh/hy_uhd_staggeredDivb
!!
!! NAME
!!
!!  hy_uhd_staggeredDivb
!!
!! SYNOPSIS
!!
!!  hy_uhd_staggeredDivb( integer(IN) :: blockID,
!!                        real(IN)    :: dt,
!!                        real(IN)    :: del(MDIM),
!!                        integer(IN) :: blkLimits,
!!                        integer(IN) :: blkLimitsGC,
!!                        logical(IN) :: halfTimeAdvance)
!!
!!
!! ARGUMENTS
!!
!!   blockID      - current block ID
!!   dt           - timestep
!!   del          - deltas in {x,y,z} directions
!!   blkLimits    - an array that holds the lower and upper indices of the section
!!                  of block without the guard cells
!!   blkLimitsGC  - an array that holds the lower and upper indices of the section
!!                  of block with the guard cells
!!   halfTimeAdvance - a flag to determine if n+1/2 time stepping is to be performed
!!
!!
!! DESCRIPTION
!!
!!   This routine computes divergence-free magnetic fields on a staggered grid.
!!   It uses the dualism relationship (Balsara & Spicer) between electric fields and the
!!   high-order Godunov fluxes and solves a set of induction equations in a discrete
!!   sense. Such calculated magnetic fields satisfy the divergence-free constraint
!!   at cell faces. The cell centered field values can be obtained by taking averages 
!!   of the cell face values.
!!
!!
!! REFERENCES
!!
!!  * Balsara & Spicer, JCP, 149:270-292, 1999
!!  * Lee, D., Ph.D. Dissertation, Univ. of MD, 2006
!!  * Lee, D. and Deane, A., "An Unsplit Staggered Mesh Scheme for Multidimensional
!!                            Magnetohydrodynamics", 228 (2009), 952-975, JCP
!!
!!*** 

!!REORDER(4): B[xyz], E, U

Subroutine hy_uhd_staggeredDivb(blockID,dt,del,blkLimits,blkLimitsGC,halfTimeAdvance)

  use Grid_interface,    ONLY : Grid_releaseBlkPtr,Grid_getBlkPtr, Grid_getCellCoords
  use Hydro_data,        ONLY : hy_energyFixSwitch, hy_geometry

  implicit none

#include "constants.h"
#include "Flash.h"
#include "UHD.h"

  !! Arguments list ----------------------------------
  integer, intent(IN) :: blockID
  real,    intent(IN) :: dt
  real,    dimension(MDIM),   intent(IN) :: del
  integer, dimension(LOW:HIGH,MDIM), intent(IN) :: blkLimits,blkLimitsGC
  logical, intent(IN) :: halfTimeAdvance
  !! -------------------------------------------------

#if NFACE_VARS > 0
#if NDIM > 1

  integer :: i,j,k,imin,imax,iskip,jmin,jmax,jskip,kmin,kmax
  real    :: dx, dy, dz
  real    :: facex,facey,facez, Ap, Am, dV, rc, rp, rm
  real, pointer, dimension(:,:,:,:) :: E,Bx,By,Bz

#ifdef FIXEDBLOCKSIZE  
  real, dimension(GRID_IHI_GC) :: xCenter
#else  
  real, dimension(blkLimitsGC(HIGH,IAXIS)) :: xCenter
#endif   



  !! Set deltas
  dx = del(DIR_X)
  dy = del(DIR_Y)
  dz = del(DIR_Z)

  imin = blkLimits(LOW, IAXIS)
  imax = blkLimits(HIGH,IAXIS)+1
  jmin = blkLimits(LOW, JAXIS)
  jmax = blkLimits(HIGH,JAXIS)+1
  if (NDIM == 2) then
     kmin = 1
     kmax = 1
  elseif (NDIM == 3) then
     kmin = blkLimits(LOW, KAXIS)
     kmax = blkLimits(HIGH,KAXIS)+1
  endif


  !! Get block pointers
  call Grid_getBlkPtr(blockID,E, SCRATCH)
  call Grid_getBlkPtr(blockID,Bx,FACEX)
  call Grid_getBlkPtr(blockID,By,FACEY)
  if (NDIM == 3) call Grid_getBlkPtr(blockID,Bz,FACEZ)


  if (hy_geometry /= CARTESIAN) then
     !get coord info will use this for Areas and Volumes
     call Grid_getCellCoords(IAXIS,blockID, CENTER,.true.,xCenter, blkLimitsGC(HIGH,IAXIS))
  endif

  if (halfTimeAdvance) then
     Bx(MAGI_FACE_VAR,:,:,:) = 0.
     By(MAGI_FACE_VAR,:,:,:) = 0.
     if (NDIM == 3) Bz(MAGI_FACE_VAR,:,:,:) = 0.
  endif

  ! Set geometric factors to be unity for Cartesian;
  ! They will be set differently otherwise.
  if (hy_geometry == CARTESIAN) then
     rc = 1.
     rp = 1.
     rm = 1.
     Ap = 1.
     Am = 1.
     dV = dx
  endif


  !! Update begins
  do k=kmin,kmax
     do j=jmin,jmax
        do i=imin,imax
#if NDIM == 2
           if (hy_geometry == CYLINDRICAL) then
              rc = xCenter(i)            
              rp = xCenter(i) + 0.5*dx
              rm = xCenter(i) - 0.5*dx
              Ap  = abs(rp)
              Am  = abs(rm)
              dV  = abs(rc)*dx
           endif

           if (j < jmax) then
              facex = Bx(MAG_FACE_VAR,i,j,k)&
                   -dt/dy*( E(EZ_SCRATCH_GRID_VAR,i,j+1,k) &
                           -E(EZ_SCRATCH_GRID_VAR,i,j,  k))
              if (halfTimeAdvance) then
                 Bx(MAGI_FACE_VAR,i,j,k) = facex
              else
                 Bx(MAG_FACE_VAR,i,j,k) = facex
              endif
           endif
           if (i < imax) then
              facey = By(MAG_FACE_VAR,i,j,k)&
                   -dt/dV*(-Ap*E(EZ_SCRATCH_GRID_VAR,i+1,j,k) &
                           +Am*E(EZ_SCRATCH_GRID_VAR,i,  j,k))
              if (halfTimeAdvance) then
                 By(MAGI_FACE_VAR,i,j,k) = facey
              else
                 By(MAG_FACE_VAR,i,j,k) = facey
              endif
           endif
#elif NDIM == 3
           if ((j < jmax) .and. (k < kmax)) then
              facex = Bx(MAG_FACE_VAR,i,j,k)&
                  -dt/dy*( E(EZ_SCRATCH_GRID_VAR,i,j+1,k  ) &
                          -E(EZ_SCRATCH_GRID_VAR,i,j,  k  ))&
                  -dt/dz*(-E(EY_SCRATCH_GRID_VAR,i,j,  k+1) &
                          +E(EY_SCRATCH_GRID_VAR,i,j,  k  ))
              if (halfTimeAdvance) then
                 Bx(MAGI_FACE_VAR,i,j,k) = facex
              else
                 Bx(MAG_FACE_VAR,i,j,k) = facex
              endif
           endif

           if ((i < imax) .and. (k < kmax)) then
              facey = By(MAG_FACE_VAR,i,j,k)&
                  -dt/dx*(-E(EZ_SCRATCH_GRID_VAR,i+1,j,k  ) &
                          +E(EZ_SCRATCH_GRID_VAR,i,  j,k  ))&
                  -dt/dz*( E(EX_SCRATCH_GRID_VAR,i,  j,k+1) &
                          -E(EX_SCRATCH_GRID_VAR,i,  j,k  ))
              if (halfTimeAdvance) then
                 By(MAGI_FACE_VAR,i,j,k) = facey
              else
                 By(MAG_FACE_VAR,i,j,k) = facey
              endif
           endif

           if ((i < imax) .and. (j < jmax)) then
              facez = Bz(MAG_FACE_VAR,i,j,k)&
                  -dt/dy*(-E(EX_SCRATCH_GRID_VAR,i,  j+1,k) &
                          +E(EX_SCRATCH_GRID_VAR,i,  j,  k))&
                  -dt/dx*( E(EY_SCRATCH_GRID_VAR,i+1,j,  k) &
                          -E(EY_SCRATCH_GRID_VAR,i,  j,  k))
              if (halfTimeAdvance) then
                 Bz(MAGI_FACE_VAR,i,j,k) = facez
              else
                 Bz(MAG_FACE_VAR,i,j,k) = facez
              endif
           endif
#endif
        enddo
     enddo
  enddo

  call Grid_releaseBlkPtr(blockID,E, SCRATCH)
  call Grid_releaseBlkPtr(blockID,Bx,FACEX)
  call Grid_releaseBlkPtr(blockID,By,FACEY)
  if (NDIM == 3) call Grid_releaseBlkPtr(blockID,Bz,FACEZ)

#endif
!end of #if NDIM > 1
#endif
!end of #if NFACE_VARS > 0
End Subroutine hy_uhd_staggeredDivb
