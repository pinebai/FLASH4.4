!!****if* source/physics/Hydro/HydroMain/split/PPM/Hydro_detectShock
!!
!!
!! NAME
!!
!!  Hydro_detectShock
!!
!!
!! SYNOPSIS
!!  call Hydro_detectShock(real(IN),pointer :: solnData(:,:,:,:), 
!!                     real(OUT)            :: shock(:,:,:), 
!!                     integer(IN)          :: blkLimits(2,MDIM),
!!                     integer(IN)          :: blkLimitsGC(2,MDIM),
!!                     integer(IN)          :: guardCells(MDIM),
!!                     real(IN)             :: primaryCoord(blkLimitsGC(LOW,IAXIS):blkLimitsGC(HIGH,IAXIS)),
!!                     real(IN)             :: secondCoord(blkLimitsGC(LOW,IAXIS):blkLimitsGC(HIGH,IAXIS)),
!!                     real(IN)             :: thirdCoord(blkLimitsGC(LOW,IAXIS):blkLimitsGC(HIGH,IAXIS)))
!!
!! DESCRIPTION
!!
!!  Hydro_detectShock is a multidimensional shock detection algorithm.  This
!!  is currently used by the burner to cut off burning in shocks if desired.
!!  This is helpful for detonations in some circumstances.
!!
!!
!! ARGUMENTS
!!
!!  solnData --     pointer to a block of solution data to inspect
!!  shock --        an array indicating if there is a shock in a zone (=1) =0 otherwise
!!  blkLimits  --   index limits of block, interior cells only; see Grid_getBlkIndexLimits 
!!  blkLimitsGC --  index limits of block, including guard cells.   
!!  guardCells  --  number of layers of guard cells which output is requested in addition
!!                  to interior cells
!!  primaryCoord -- x coordinate of solnData, i.e., coordinate in IAXIS direction, including guardcells
!!  secondCoord --  y coordinate of solnData, i.e., coordinate in JAXIS direction, including guardcells
!!  thirdCoord --   z coordinate of solnData, i.e., coordinate in KAXIS direction, including guardcells
!!
!!  NOTES
!!
!!  The guard cells need to be filled before calling this routine.
!!  The pressure needs to be updated before calling this routine.
!!  
!!  SEE ALSO
!!
!!  Grid_getBlkIndexLimits
!!
!!***

! solnData depends on the ordering on unk
!!REORDER(4): solnData



#ifdef DEBUG_ALL
#define DEBUG_HYDRO
#endif

subroutine Hydro_detectShock(solnData, shock, blkLimits, blkLimitsGC, &
                             guardCells, &
                             primaryCoord,secondCoord,thirdCoord)

  use Hydro_data, ONLY :  hy_smallu, hy_geometry, hy_dp_sh_md
  use Driver_interface, ONLY : Driver_abortFlash

  implicit none
#include "constants.h"
#include "Flash.h"

  integer, intent(IN), dimension(2,MDIM) :: blkLimits, blkLimitsGC
  integer, intent(IN) :: guardCells(MDIM)
  real, pointer :: solnData(:,:,:,:) 
#ifdef FIXEDBLOCKSIZE
  real, intent(out),dimension(GRID_ILO_GC:GRID_IHI_GC,&
                              GRID_JLO_GC:GRID_JHI_GC,&
                              GRID_KLO_GC:GRID_KHI_GC):: shock
  real,intent(IN),dimension(GRID_ILO_GC:GRID_IHI_GC) :: primaryCoord
  real,intent(IN),dimension(GRID_JLO_GC:GRID_JHI_GC) :: secondCoord
  real,intent(IN),dimension(GRID_KLO_GC:GRID_KHI_GC) :: thirdCoord
  real, dimension(MDIM,GRID_ILO_GC:GRID_IHI_GC,&
                  GRID_JLO_GC:GRID_JHI_GC,&
                  GRID_KLO_GC:GRID_KHI_GC):: div_v
#else
  real,intent(out), dimension(blkLimitsGC(LOW,IAXIS):blkLimitsGC(HIGH,IAXIS),&
                               blkLimitsGC(LOW,JAXIS):blkLimitsGC(HIGH,JAXIS),&
                               blkLimitsGC(LOW,KAXIS):blkLimitsGC(HIGH,KAXIS)) :: shock
  real,intent(IN),dimension(blkLimitsGC(LOW,IAXIS):blkLimitsGC(HIGH,IAXIS)) :: primaryCoord
  real,intent(IN),dimension(blkLimitsGC(LOW,JAXIS):blkLimitsGC(HIGH,JAXIS)) :: secondCoord
  real,intent(IN),dimension(blkLimitsGC(LOW,KAXIS):blkLimitsGC(HIGH,KAXIS)) :: thirdCoord
  real, dimension(MDIM,blkLimitsGC(LOW,IAXIS):blkLimitsGC(HIGH,IAXIS),&
                  blkLimitsGC(LOW,JAXIS):blkLimitsGC(HIGH,JAXIS),&
                  blkLimitsGC(LOW,KAXIS):blkLimitsGC(HIGH,KAXIS)) :: div_v
#endif

  real :: smlusq

  integer :: i, j, k
  integer :: ib, ie, jb, je, kb, ke

  real :: ddpx, p0x, p1x
  real :: ddxux, ddxuxm
  
  real :: ddpy, p0y, p1y
  real :: ddyuy, ddyuym
      
  real :: ddpz, p0z, p1z
  real :: ddzuz, ddzuzm
  
  real :: factrx, factry, factrz
  real :: deenom, ppre, ppst

  real :: flag
  real :: divu, shockd

!------------------------------------------------------------------------------

#ifdef DEBUG_HYDRO
  if (guardCells .GE. NGUARD) then
     call Driver_abortFlash("Hydro_detectShock: too many guard cell layers requested")
  end if
#endif

! get a pointer to all the data for the current block

  ib = max(blkLimitsGC(LOW,IAXIS)+1,    blkLimits(LOW,IAXIS) -     guardCells(IAXIS) )
  ie = min(blkLimitsGC(HIGH,IAXIS)-1,   blkLimits(HIGH,IAXIS)+     guardCells(IAXIS) )
  jb = max(blkLimitsGC(LOW,JAXIS)+K2D,  blkLimits(LOW,JAXIS) - K2D*guardCells(JAXIS) )
  je = min(blkLimitsGC(HIGH,JAXIS)-K2D, blkLimits(HIGH,JAXIS)+ K2D*guardCells(JAXIS) )
  kb = max(blkLimitsGC(LOW,KAXIS)+K3D,  blkLimits(LOW,KAXIS) - K3D*guardCells(KAXIS) )
  ke = min(blkLimitsGC(HIGH,KAXIS)-K3D, blkLimits(HIGH,KAXIS)+ K3D*guardCells(KAXIS) )


  smlusq = hy_smallu**2


#ifdef DEBUG_HYDRO
  if((hy_geometry /= CARTESIAN).and.&
     (hy_geometry /= CYLINDRICAL).and.&
     (hy_geometry /= SPHERICAL).and.&
     (hy_geometry /= POLAR))call Driver_abortFlash("shockDetect:wrong geometry")
  if((NDIM ==1).and.(hy_geometry/=CARTESIAN).and.(hy_geometry/=SPHERICAL))&
       call Driver_abortFlash("SHOCK DETECT : 1d cylindrical not supported")
  if((NDIM ==3).and.(hy_geometry==POLAR))&
       call Driver_abortFlash("SHOCK DETECT : 3d polar not supported")
#endif

#if NDIM == 1
! compute the velocity divergence
  if (hy_geometry == CARTESIAN) then
     do i = ib, ie
        
        div_v(IAXIS,i,1,1) = (solnData(VELX_VAR,i+1,1,1) -    &
                             solnData(VELX_VAR,i-1,1,1)) /   &
                             (primaryCoord(i+1) - primaryCoord(i-1))
        
     enddo
        
  else if (hy_geometry == SPHERICAL) then
        
     do i = ib, ie
           
        ! divergence is (1/r**2) (d/dr) (r**2 v_r)
        
        div_v(IAXIS,i,1,1) =                 &
             (primaryCoord(i+1)**2*solnData(VELX_VAR,i+1,1,1) - &
             primaryCoord(i-1)**2*solnData(VELX_VAR,i-1,1,1)) / &
             (primaryCoord(i)**2*(primaryCoord(i+1) - primaryCoord(i-1)))
           
     enddo
        
  endif
  
#elif NDIM == 2
  if (hy_geometry == CARTESIAN) then
     do j = jb, je
        do i = ib, ie
              
           div_v(IAXIS,i,j,1) =                 &
                (solnData(VELX_VAR,i+1,j,1) -    & 
                solnData(VELX_VAR,i-1,j,1)) /   &
                (primaryCoord(i+1) - primaryCoord(i-1))
              
           div_v(JAXIS,i,j,1) =                 &
                (solnData(VELY_VAR,i,j+1,1) -    &
                solnData(VELY_VAR,i,j-1,1)) /   &
                (secondCoord(j+1) - secondCoord(j-1))
              
        enddo
     enddo
        
  elseif (hy_geometry == CYLINDRICAL) then
        
        ! 2-d cylindrical geometry (r,z)
        
     do j = jb, je
        do i = ib, ie
              
           ! divergence is (1/r) (d/dr) (v_r) + (d/dz) (v_z)
           
           div_v(IAXIS,i,j,1) =                 &
                (primaryCoord(i+1)*solnData(VELX_VAR,i+1,j,1) -    &
                primaryCoord(i-1)*solnData(VELX_VAR,i-1,j,1)) /   &
                (primaryCoord(i)*(primaryCoord(i+1) - primaryCoord(i-1)))
              
           div_v(JAXIS,i,j,1) =                 &
                (solnData(VELY_VAR,i,j+1,1) -    & 
                solnData(VELY_VAR,i,j-1,1)) /   &
                (secondCoord(j+1) - secondCoord(j-1))  
           
        enddo
     enddo
        
  elseif (hy_geometry == POLAR) then
        
     ! 2-d polar geometry (r,phi)
        
     do j = jb, je
        do i = ib, ie
              
           ! divergence is (1/r) (d/dr) (v_r) + (1/r) (d/df) (v_f)
           
           div_v(IAXIS,i,j,1) =                 &
                (primaryCoord(i+1)*solnData(VELX_VAR,i+1,j,1) -    &
                primaryCoord(i-1)*solnData(VELX_VAR,i-1,j,1)) /   &
                (primaryCoord(i)*(primaryCoord(i+1) - primaryCoord(i-1)))
           
           div_v(JAXIS,i,j,1) =                 &
                (solnData(VELY_VAR,i,j+1,1) -    & 
                solnData(VELY_VAR,i,j-1,1)) /   &
                (primaryCoord(i)*(secondCoord(j+1) - secondCoord(j-1)))
           
        enddo
     enddo
        
  elseif (hy_geometry == SPHERICAL) then
        
     ! spherical geometry (r,theta)
     
     do j = jb, je
        do i = ib, ie
           
           ! divergence is (1/r**2) (d/dr) (r**2 v_r) 
           ! + (1/(r sin(t)) ) (d/dt)(sin(t) v_t)
           
           div_v(IAXIS,i,j,1) = &
                (primaryCoord(i+1)*primaryCoord(i+1)&
                *solnData(VELX_VAR,i+1,j,1) - &
                primaryCoord(i-1)*primaryCoord(i-1)&
                *solnData(VELX_VAR,i-1,j,1)) / &
                (primaryCoord(i)*primaryCoord(i)*(primaryCoord(i+1) &
                - primaryCoord(i-1)))
           
           div_v(JAXIS,i,j,1) = &
                (sin(secondCoord(j+1))*solnData(VELY_VAR,i,j+1,1) - &
                sin(secondCoord(j-1))*solnData(VELY_VAR,i,j-1,1)) / &
                (primaryCoord(i)*sin(secondCoord(j))*&
               (secondCoord(j+1) - secondCoord(j-1)))
              
        enddo
     enddo
     
  endif
  
#else
  
  if (hy_geometry == CARTESIAN) then 
     
     do k = kb,ke
        do j = jb,je
           do i = ib,ie
              
              div_v(IAXIS,i,j,k) =                 &
                   (solnData(VELX_VAR,i+1,j,k) -    &
                   solnData(VELX_VAR,i-1,j,k)) /   &
                   (primaryCoord(i+1) - primaryCoord(i-1))           
              
              div_v(JAXIS,i,j,k) =                 &
                   (solnData(VELY_VAR,i,j+1,k) -    &
                   solnData(VELY_VAR,i,j-1,k)) /   &
                   (secondCoord(j+1) - secondCoord(j-1)) 
              
              div_v(KAXIS,i,j,k) =                 &
                   (solnData(VELZ_VAR,i,j,k+1) -    &
                   solnData(VELZ_VAR,i,j,k-1)) /   &
                   (thirdCoord(k+1) - thirdCoord(k-1))
              
           enddo
        enddo
     enddo
     
  else if (hy_geometry == CYLINDRICAL) then
     
     ! 3-d cylindrical geometry (r,z,phi)
     
     do k = kb,ke
        do j = jb,je
           do i = ib,ie
              
              ! divergence is (1/r) (d/dr) (v_r) 
              ! + (d/dz) (v_z) + (1/r) (d/df) (v_f)
              
              div_v(IAXIS,i,j,k) =                        &
                   (primaryCoord(i+1)*solnData(VELX_VAR,i+1,j,k) -    &
                   primaryCoord(i-1)*solnData(VELX_VAR,i-1,j,k)) /   &
                   (primaryCoord(i)*(primaryCoord(i+1) - primaryCoord(i-1)))
              
              div_v(JAXIS,i,j,k) =                 &
                   (solnData(VELY_VAR,i,j+1,k) -    &
                   solnData(VELY_VAR,i,j-1,k)) /   &
                   (secondCoord(j+1) - secondCoord(j-1))
              
              div_v(KAXIS,i,j,k) =                 &
                   (solnData(VELZ_VAR,i,j,k+1) -    &
                   solnData(VELZ_VAR,i,j,k-1)) /   &
                   (primaryCoord(i)*(thirdCoord(k+1) - thirdCoord(k-1)))
              
           enddo
        enddo
     enddo

  else if (hy_geometry == SPHERICAL) then
     
     ! 3-d spherical geometry (r,theta,phi)
     
     do k = kb,ke
        do j = jb,je
           do i = ib,ie
              
              ! divergence is (1/r**2) (d/dr) (r**2 v_r) 
              ! + (1/(r sin(t)) ) (d/dt)(sin(t) v_t)
              ! + (1/(r sin(t)) ) (d/df) (v_f)
              
              div_v(IAXIS,i,j,k) = &
                   (primaryCoord(i+1)*primaryCoord(i+1)&
                   *solnData(VELX_VAR,i+1,j,k) - &
                   primaryCoord(i-1)*primaryCoord(i-1)&
                   *solnData(VELX_VAR,i-1,j,k)) / &
                   (primaryCoord(i)*primaryCoord(i)*(primaryCoord(i+1) &
                   - primaryCoord(i-1)))
              
              div_v(JAXIS,i,j,k) = &
                (sin(secondCoord(j+1))*solnData(VELY_VAR,i,j+1,k) - &
                sin(secondCoord(j-1))*solnData(VELY_VAR,i,j-1,k)) / &
                (primaryCoord(i)*sin(secondCoord(j))*&
               (secondCoord(j+1) - secondCoord(j-1)))
              
              div_v(KAXIS,i,j,k) =                 &
                   (solnData(VELZ_VAR,i,j,k+1) -    &
                   solnData(VELZ_VAR,i,j,k-1)) /   &
                   (primaryCoord(i)*sin(secondCoord(j))*(thirdCoord(k+1) - thirdCoord(k-1)))
              
           enddo
        enddo
     enddo
  endif
  
#endif
  
  ! shock detection and direction of shock propagation
  
  do k = kb,ke
     do j = jb,je
        do i = ib,ie
           
! interface centred jumps along given dimension needed for calculation 
! of direction of shock propagation

! first dimension
           ddpx = solnData(PRES_VAR,i+1,j,k) - solnData(PRES_VAR,i-1,j,k)

           if (ddpx .LT. 0.e0) then
              p0x  = solnData(PRES_VAR,i+1,j,k)
              p1x  = solnData(PRES_VAR,i-1,j,k)
           else
              p0x  = solnData(PRES_VAR,i-1,j,k)
              p1x  = solnData(PRES_VAR,i+1,j,k)
           end if

! look for compression, if vx(i+1) - vx(i-1) < 0, then we are compressing in x
           ddxux = solnData(VELX_VAR,i+1,j,k) - solnData(VELX_VAR,i-1,j,k)
           ddxuxm = min(0.e0, ddxux )

! second dimension
           if(NDIM .GT. 1) then
             ddpy = solnData(PRES_VAR,i,j+1,k) - solnData(PRES_VAR,i,j-1,k)

             if (ddpy .LT. 0.e0) then
                p0y  = solnData(PRES_VAR,i,j+1,k)
                p1y  = solnData(PRES_VAR,i,j-1,k)     
             else
                p0y  = solnData(PRES_VAR,i,j-1,k)
                p1y  = solnData(PRES_VAR,i,j+1,k)
             end if
               
! look for compression, if vy(i+1) - vy(i-1) < 0, then we are compressing in y
             ddyuy  = solnData(VELY_VAR,i,j+1,k) - solnData(VELY_VAR,i,j-1,k)
             ddyuym = min(0.e0, ddyuy )
          else
             p0y    = 0.e0
             p1y    = 0.e0
             ddyuym = 0.e0
          end if
          ! third dimension
          if (NDIM .EQ. 3) then 
             
             ddpz = solnData(PRES_VAR,i,j,k+1) - solnData(PRES_VAR,i,j,k-1)
             
             if (ddpz .LT. 0.e0) then
                p0z  = solnData(PRES_VAR,i,j,k+1)
                p1z  = solnData(PRES_VAR,i,j,k-1)
             else
                p0z  = solnData(PRES_VAR,i,j,k-1)
                p1z  = solnData(PRES_VAR,i,j,k+1)
             end if
             
             ! look for compression, if vz(i+1) - vz(i-1) < 0, then we are compressing in z
             ddzuz  = solnData(VELZ_VAR,i,j,k+1) - solnData(VELZ_VAR,i,j,k-1)
             ddzuzm = min(0.e0, ddzuz )
          else
             p0z    = 0.e0
             p1z    = 0.e0
             ddzuzm = 0.e0
          end if


! We will use the factors factr(xyz) to blend data from the longitudinal 
! and transverse directions in order to estimate shock jumps, orientation, 
! speeds, etc., below.  This data blending technique generates more algebra
! here, but it eliminates any tendency to flip-flop near a 45 degree shock 
! orientation. These blending factors are based upon velocity jumps in 
! order to sense the direction of any compression, even if that compression
! is not momentarily associated with a pressure jump (as in the case of a 
! shock reflecting from a rigid wall).

! blend the pre- and post-shock values
          factrx = ddxuxm * ddxuxm
          factry = ddyuym * ddyuym
          factrz = ddzuzm * ddzuzm
          
          deenom = 1.e0 / (factrx + factry + factrz + smlusq)
          
          factry = factry * deenom
          factrz = factrz * deenom
          factrx = 1.e0 - (factry + factrz)
          
          ! compute the pre-shock and post-shock pressures across the shock front
          ppre = factrx*p0x + factry*p0y + factrz*p0z
          ppst = factrx*p1x + factry*p1y + factrz*p1z
          
          if (ppre.NE.0.0) then
             flag = hy_dp_sh_md - (ppst - ppre) / ppre
          else if (ppst==0.0) then
             flag = 0.0
          else
             flag = -1.0
          end if
          
! construct the divergence of the velocity
          if (NDIM .EQ. 1) then
             divu = div_v(IAXIS,i,j,k)
          elseif (NDIM .EQ. 2) then
             divu = div_v(IAXIS,i,j,k) + div_v(JAXIS,i,j,k)
          else
             divu = div_v(IAXIS,i,j,k) + div_v(JAXIS,i,j,k) + div_v(KAXIS,i,j,k)
          endif

! decide is we have a shock
          shockd = max(divu, flag)
          
          if (shockd .LT. 0.e0) then
             shockd = 1.e0
          else
             shockd = 0.e0
          end if
          
          shock(i,j,k) = shockd
          
       enddo
    enddo
 enddo
 
 
end subroutine Hydro_detectShock




