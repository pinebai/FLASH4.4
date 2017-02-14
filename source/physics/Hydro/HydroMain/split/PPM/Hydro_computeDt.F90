!!****if* source/physics/Hydro/HydroMain/split/PPM/Hydro_computeDt
!!
!! NAME
!!
!!  Hydro_computeDt
!!
!!
!! SYNOPSIS
!!
!!  Hydro_computeDt(integer(IN):: blockID, 
!!                  real(IN) :: x(:), 
!!                  real(IN) :: dx(:), 
!!                  real(IN) :: uxgrid(:),
!!                  real(IN) :: y(:), 
!!                  real(IN) :: dy(:), 
!!                  real(IN) :: uygrid(:), 
!!                  real(IN) :: z(:), 
!!                  real(IN) :: dz(:), 
!!                  real(IN) :: uzgrid(:), 
!!                  integer(IN) :: blkLimits(2,MDIM)
!!                  integer(IN) :: blkLimitsGC(2,MDIM)
!!                  real,pointer ::  solnData(:,:,:,:),   
!!                  real,(INOUT) ::   dtCheck, 
!!                  integer(INOUT) :: dtMinLoc(:),
!!                  real(INOUT), optional :: extraInfo)
!!
!! DESCRIPTION
!!
!!  Computes the timestep limiter for the hydrodynamical solver.  For pure
!!  hydrodynamics, the Courant-Fredrichs-Lewy criterion is used.  The sound
!!  speed is computed and together with the velocities, is used to constrain
!!  the timestep such that no information can propagate more than one zone
!!  per timestep.
!!
!!
!! ARGUMENTS
!!
!!  blockID -       local block ID
!!  x, y, z -       coordinates
!!  dx, dy, dz -    deltas in each {x, y z} directions
!!  uxgrid, uygrid, uzgrid - velocity of grid expansion in {x, y z} directions
!!  blkLimits -    the indices for the interior endpoints of the block
!!  blkLimitsGC - the indices for endpoints including the guardcells
!!  solnData -      the physical, solution data from grid
!!  dtCheck -      variable to hold timestep constraint
!!  dtMinLoc(5) -  array to hold location of cell responsible for minimum dt:
!!                 dtMinLoc(1) = i index
!!                 dtMinLoc(2) = j index
!!                 dtMinLoc(3) = k index
!!                 dtMinLoc(4) = blockID
!!                 dtMinLoc(5) = hy_meshMe
!! extraInfo    -  Driver_computeDt can provide extra info to the caller
!!                 using this argument.
!!
!!***

! solnData depends on the ordering on unk
!!REORDER(4): solnData


subroutine Hydro_computeDt (blockID,  &
                           x, dx, uxgrid, &
                           y, dy, uygrid, &
                           z, dz, uzgrid, &
                           blkLimits,blkLimitsGC,        &
                           solnData,   &
                           dtCheck, dtMinLoc,&
                           extraInfo )
     
  
#include "Flash.h"
#include "constants.h"

  use Hydro_data, ONLY : hy_useHydro, hy_geometry, hy_cfl, hy_meshMe
  use Driver_interface, ONLY : Driver_abortFlash

  implicit none


  integer, intent(IN) :: blockID 
  integer, intent(IN),dimension(LOW:HIGH,MDIM)::blkLimits,blkLimitsGC
  real,INTENT(INOUT)    :: dtCheck
  integer,INTENT(INOUT)    :: dtMinLoc(5)
  real, pointer :: solnData(:,:,:,:) 
#ifdef FIXEDBLOCKSIZE
  real, dimension(GRID_ILO_GC:GRID_IHI_GC), intent(IN) :: x, dx, uxgrid
  real, dimension(GRID_JLO_GC:GRID_JHI_GC), intent(IN) :: y, dy, uygrid
  real, dimension(GRID_KLO_GC:GRID_KHI_GC), intent(IN) :: z, dz, uzgrid
#else
  real, dimension(blkLimitsGC(LOW,IAXIS):blkLimitsGC(HIGH,IAXIS)), intent(IN) :: x, dx, uxgrid
  real, dimension(blkLimitsGC(LOW,JAXIS):blkLimitsGC(HIGH,JAXIS)), intent(IN) :: y, dy, uygrid
  real, dimension(blkLimitsGC(LOW,KAXIS):blkLimitsGC(HIGH,KAXIS)), intent(IN) :: z, dz, uzgrid
#endif
  real,OPTIONAL,intent(INOUT) :: extraInfo
  
  integer :: i, j, k, temploc(5)
  real    :: sndspd, sndspd2, delxinv, delyinv, delzinv, dt_temp, dt_ltemp

  
!==============================================================================


  dt_temp    = 0.
  temploc(:) = 0
  

  if (.NOT. hy_useHydro) return

  if (hy_geometry == CYLINDRICAL .AND. NDIM == 3) then
     if (hy_meshMe .EQ. MASTER_PE) print *, 'ERROR: 3-d cylindrical not yet supported in Hydro_computeDt'
     call Driver_abortFlash("ERROR: 3-d cylindrical not yet supported in Hydro_computeDt")
  endif

  !--------------------------------------------------------------------------
  ! 1-dimensional
  !--------------------------------------------------------------------------

#ifndef FIXEDBLOCKSIZE
  i=blkLimitsGC(HIGH,IAXIS)-blkLimitsGC(LOW,IAXIS)+1
  j=blkLimitsGC(HIGH,JAXIS)-blkLimitsGC(LOW,JAXIS)+1
  k=blkLimitsGC(HIGH,KAXIS)-blkLimitsGC(LOW,KAXIS)+1
#endif  

  if (NDIM == 1) then
     
     ! here we assume that the grid is uniformly spaced within each block
     delxinv = 1.0/dx(blkLimits(LOW,IAXIS))
     
     do i = blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS)
        
        sndspd = sqrt(solnData(GAMC_VAR,i,1,1) *   &
             solnData(PRES_VAR,i,1,1) /          & 
             solnData(DENS_VAR,i,1,1) )   
           
        dt_ltemp = (abs(solnData(VELX_VAR,i,1,1) - &
             uxgrid(i)) + sndspd) * delxinv
        if (dt_ltemp > dt_temp) then
           dt_temp    = dt_ltemp
           temploc(1) = i
           temploc(2) = 1
           temploc(3) = 1
           temploc(4) = blockID
           temploc(5) = hy_meshMe
        endif

     enddo

     !--------------------------------------------------------------------------
     ! 2-dimensional
     !--------------------------------------------------------------------------
     
     ! angular coordinates may be present

  elseif (NDIM == 2) then
     
     if (hy_geometry == CARTESIAN .OR. hy_geometry == CYLINDRICAL) then        
        
        ! the 'y' coordinate is not angular in Cartesian and cylindrical coords
        
        ! here we assume that the grid is uniformly spaced within each block
        delxinv = 1.0/dx(blkLimits(LOW,IAXIS))
        delyinv = 1.0/dy(blkLimits(LOW,JAXIS))
        do j = blkLimits(LOW,JAXIS), blkLimits(HIGH,JAXIS)
           do i = blkLimits(LOW,IAXIS), blkLimits(HIGH,IAXIS)

              sndspd2 = solnData(GAMC_VAR,i,j,1) * &
                   solnData(PRES_VAR,i,j,1) /       &
                   solnData(DENS_VAR,i,j,1) 
                 
              !!              dt_ltemp = abs(solnData(VELX_VAR,i,j,1) - &
              !!                             uxgrid(i)) * delxinv + &
              !!                         abs(solnData(VELY_VAR,i,j,1) - &
              !!                             uygrid(j)) * delyinv + &
              !!                         sqrt(sndspd2 * (delxinv**2 + delyinv**2))
              
              !!New time step calculation.
              !!
              !!dt = C min (dx/(|vx|+c), dy/(|vy|+c), dz/(|vz|+c))
              !!
              !!In the previous form, for the limiting case of vx=vy=vz=c and dx=dy=dz,
              !!the timestep was smaller than needed by 1.7 in 2d, and by 2.4 in 3d.
              !!
                 
              dt_ltemp = max( (abs(solnData(VELX_VAR,i,j,1) - uxgrid(i)) + &
                   sqrt(sndspd2)) * delxinv,   &
                   (abs(solnData(VELY_VAR,i,j,1) - uygrid(j)) + &
                   sqrt(sndspd2)) * delyinv )
                 
              if (dt_ltemp > dt_temp) then
                 dt_temp    = dt_ltemp
                 temploc(1) = i
                 temploc(2) = j
                 temploc(3) = 1
                 temploc(4) = blockID
                 temploc(5) = hy_meshMe
              endif

           enddo
        enddo

     else                         


        ! Angular
        
        do j = blkLimits(LOW,JAXIS), blkLimits(HIGH,JAXIS)
           do i = blkLimits(LOW,IAXIS), blkLimits(HIGH,IAXIS)
              
              ! get delxinv and delyinv for each zone -- so things could be uniformly spaced
              ! or logarithmically spaced here.
              delxinv = 1.0/dx(i)
              delyinv = 1.0/(x(i)*dy(j))
                 
              sndspd2 = solnData(GAMC_VAR,i,j,1) * &
                   solnData(PRES_VAR,i,j,1) /       &
                   solnData(DENS_VAR,i,j,1) 
                 
              !!              dt_ltemp = abs(solnData(VELX_VAR,i,j,1) -               &
              !!                             uxgrid(i)) * delxinv +       &
              !!                         abs(solnData(VELY_VAR,i,j,1) -             &
              !!                             uygrid(j)) * delyinv +       &
              !!                         sqrt(sndspd2 * (delxinv**2 + delyinv**2))
                 
              !!New time step calculation.
              !!
              !!dt = C min (dx/(|vx|+c), dy/(|vy|+c), dz/(|vz|+c))
              !!
              !!In the previous form, for the limiting case of vx=vy=vz=c and dx=dy=dz,
              !!the timestep was smaller than needed by 1.7 in 2d, and by 2.4 in 3d.
              !!
                 
              dt_ltemp = max( (abs(solnData(VELX_VAR,i,j,1) - uxgrid(i)) + &
                   sqrt(sndspd2)) * delxinv,   &
                   (abs(solnData(VELY_VAR,i,j,1) - uygrid(j)) + &
                   sqrt(sndspd2)) * delyinv )
                 
              if (dt_ltemp > dt_temp) then
                 dt_temp    = dt_ltemp
                 temploc(1) = i
                 temploc(2) = j
                 temploc(3) = 1
                 temploc(4) = blockID
                 temploc(5) = hy_meshMe
              endif

           enddo
        enddo

     endif



     !--------------------------------------------------------------------------
     ! 3-dimensional
     !--------------------------------------------------------------------------
     
     ! here we only allow Cartesian and spherical right now
     
  else
     
     
     if (hy_geometry == CARTESIAN) then         
        
        ! compute the inverse of the grid spacing.  Here we assume that the mesh is 
        ! uniformly spaced within a block
        delxinv = 1.0/dx(blkLimits(LOW,IAXIS))
        delyinv = 1.0/dy(blkLimits(LOW,JAXIS))       
        delzinv = 1.0/dz(blkLimits(LOW,KAXIS))
        
        do k = blkLimits(LOW,KAXIS), blkLimits(HIGH,KAXIS)        
           do j = blkLimits(LOW,JAXIS), blkLimits(HIGH,JAXIS)                         
              do i = blkLimits(LOW,IAXIS), blkLimits(HIGH,IAXIS)

                 sndspd2 =  solnData(GAMC_VAR,i,j,k) * &
                      solnData(PRES_VAR,i,j,k) /  &
                      solnData(DENS_VAR,i,j,k) 
                    
                 !! Updated/fixed time step below. This was provious, wrong one.
                 !!                 dt_ltemp = abs(solnData(VELX_VAR,i,j,k) - &
                 !!                                uxgrid(i)) * delxinv + &
                 !!                            abs(solnData(VELY_VAR,i,j,k) - &
                 !!                                uygrid(j)) * delyinv + &
                 !!                            abs(solnData(VELZ_VAR,i,j,k) - &
                 !!                                uzgrid(k)) * delzinv + &
                 !!                          sqrt(sndspd2 * (delxinv**2 + delyinv**2 + delzinv**2))
                 
                 !!New time step calculation.
                 !!
                 !!dt = C min (dx/(|vx|+c), dy/(|vy|+c), dz/(|vz|+c))
                 !!
                 !!In the previous form, for the limiting case of vx=vy=vz=c and dx=dy=dz,
                 !!the timestep was smaller than needed by 1.7 in 2d, and by 2.4 in 3d.
                 !!
                    
                    
                 dt_ltemp = max( (abs(solnData(VELX_VAR,i,j,k) - uxgrid(i)) + &
                      sqrt(sndspd2)) * delxinv,   &
                      (abs(solnData(VELY_VAR,i,j,k) - uygrid(j)) + &
                      sqrt(sndspd2)) * delyinv ,  &
                      (abs(solnData(VELZ_VAR,i,j,k) - uzgrid(k)) + &
                      sqrt(sndspd2)) * delzinv )
                    
                 if (dt_ltemp > dt_temp) then
                    dt_temp    = dt_ltemp
                    temploc(1) = i
                    temploc(2) = j
                    temploc(3) = k
                    temploc(4) = blockID
                    temploc(5) = hy_meshMe
                 endif

              enddo
           enddo
        enddo

        
     elseif (hy_geometry == SPHERICAL) then
        
        do k = blkLimits(LOW,KAXIS), blkLimits(HIGH,KAXIS)        
           do j = blkLimits(LOW,JAXIS), blkLimits(HIGH,JAXIS)           
              do i = blkLimits(LOW,IAXIS), blkLimits(HIGH,IAXIS)

                 ! allow for the grid to be non-uniformly spaced within a block
                 delzinv = 1.0/(x(i)*sin(y(j))*dz(k))
                 delyinv = 1.0/(x(i)*dy(j))
                 delxinv = 1.0/dx(i)
                    
                 sndspd2 = solnData(GAMC_VAR,i,j,k) * &
                      solnData(PRES_VAR,i,j,k) /  &
                      solnData(DENS_VAR,i,j,k) 
                    
                 !!                 dt_ltemp = abs(solnData(VELX_VAR,i,j,k) - &
                 !!                               uxgrid(i)) * delxinv + &
                 !!                            abs(solnData(VELY_VAR,i,j,k) - &
                 !!                                uygrid(j)) * delyinv + &
                 !!                            abs(solnData(VELZ_VAR,i,j,k) - &
                 !!                                uzgrid(k)) * delzinv + &
                 !!                         sqrt(sndspd2 * (delxinv**2 + delyinv**2 + delzinv**2))
                    
                 !!New time step calculation.
                 !!
                 !!dt = C min (dx/(|vx|+c), dy/(|vy|+c), dz/(|vz|+c))
                 !!
                 !!In the previous form, for the limiting case of vx=vy=vz=c and dx=dy=dz,
                 !!the timestep was smaller than needed by 1.7 in 2d, and by 2.4 in 3d.
                 !!
                    
                 dt_ltemp = max( (abs(solnData(VELX_VAR,i,j,k) - uxgrid(i)) + &
                      sqrt(sndspd2)) * delxinv,   &
                      (abs(solnData(VELY_VAR,i,j,k) - uygrid(j)) + &
                      sqrt(sndspd2)) * delyinv ,  &
                      (abs(solnData(VELZ_VAR,i,j,k) - uzgrid(k)) + &
                      sqrt(sndspd2)) * delzinv )
                 
                 if (dt_ltemp > dt_temp) then
                    dt_temp    = dt_ltemp
                    temploc(1) = i
                    temploc(2) = j
                    temploc(3) = k
                    temploc(4) = blockID
                    temploc(5) = hy_meshMe
                 endif
                 
              enddo
           enddo
        enddo

     else

        call Driver_abortFlash("ERROR: geometry not supported in Hydro_computeDt")
        
     endif
     
  endif

  if (dt_temp == 0.) dt_temp = 1.d-100

  dt_temp = hy_cfl / dt_temp
  if (dt_temp < dtCheck) then
     dtCheck = dt_temp
     dtMinLoc = temploc
  endif
  if(dtCheck <= 0.0) call Driver_abortFlash("[Hydro]: Computed dt is not positive! Aborting!")
  
  !! Set default output to be null if not needed.
  if (present(extraInfo)) extraInfo = 0.

  return
end subroutine Hydro_computeDt


