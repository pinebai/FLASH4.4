!!****if* source/Particles/ParticlesMapping/Quadratic/pt_mapFromMeshQuadratic
!!
!! NAME
!!
!!  pt_mapFromMeshQuadratic
!!
!! SYNOPSIS
!!
!!  pt_mapFromMeshQuadratic(integer, INTENT(in)    :: numAttrib,
!!                        integer, INTENT(in)    :: attrib(2,numAttrib),
!!                           real, INTENT(in)    :: pos(MDIM),
!!                           real, INTENT(in)    :: bndBox(2,MDIM),
!!                           real, INTENT(in)    :: deltaCell(MDIM),
!!                           real, pointer       :: solnVec(:,:,:,:),
!!                           real, INTENT(OUT)   :: partAttribVec(numAttrib))
!!
!! DESCRIPTION
!!
!!  Routine to map a list of  quantities defined on the mesh onto a single
!!  particle position.
!!
!!  Currently volume elements are assumed to be Cartesian or
!!  2D axisymmetric (r-z).
!!
!!  This version does quadratic interpolation onto particle
!!  positions (rather than the standard linear particle-mesh
!!  mapping).
!!  The routines assumes that the input arguement attrib contains the 
!!  list of particle attributes that need updating, and the corresponding
!!  mesh variables. The following constants are define in "Particles.h"
!!  which must be included in the routine
!!
!!     PART_DS_IND  : index into particle data structure   
!!     GRID_DS_IND  : index into the grid data structure
!!  As and example, if it is desired to map
!!  all three directional velocities, and the density then
!!  numAttrib is 4, and the array attrib will have the following values
!!     attrib(PART_DS_IND,1)=VELX_PART_PROP, attrib(GRID_DS_IND,1)=VELX_VAR
!!     attrib(PART_DS_IND,2)=VELY_PART_PROP, attrib(GRID_DS_IND,2)=VELY_VAR
!!     attrib(PART_DS_IND,3)=VELZ_PART_PROP, attrib(GRID_DS_IND,3)=VELZ_VAR
!!     attrib(PART_DS_IND,4)=DENS_PART_PROP, attrib(GRID_DS_IND,4)=DENS_VAR
!!  Here the order of the attributes is completely unimportant, only the 
!!  map between the particle property and the grid variable should be on the
!!  same index
!!
!! ARGUMENTS
!!
!!   numAttrib : number of attributes to update
!!   attrib    : list containing the attributes and the corresponding grid
!!               data structure index
!!   pos       : The physical coordinates of the particle
!!   bndBox    : the bounding box of the block containing the particle
!!   deltaCell : the dx,dy,dz of the block
!!   solnVec   : Pointer to the solution data block in which particle is found
!!   partAttribVec  : calculated particle attribute values 
!!
!! NOTES
!!
!!  The PART_DS_IND information in attrib is not used by this subroutine.
!!  It is the caller's responsibility to copy the interpolation results
!!  returned in partAttribVec into the appropriate particle property
!!  slots in accordance with attrib(PART_DS_IND,:).
!!
!!***

!!REORDER(4):solnVec

subroutine pt_mapFromMeshQuadratic (numAttrib, attrib, pos, bndBox,&
     deltaCell,solnVec, partAttribVec)
  
  use Particles_data, ONLY : pt_geometry, pt_str_geometry
  use Driver_interface, ONLY : Driver_abortFlash

  implicit none

#include "constants.h"
#include "Flash.h"
#include "Particles.h"

  integer, INTENT(in) :: numAttrib
  integer, dimension(2, numAttrib),intent(IN) :: attrib
  real,dimension(MDIM), INTENT(in)    :: pos,deltaCell
  real, dimension(LOW:HIGH,MDIM), intent(IN) :: bndBox
  real, pointer       :: solnVec(:,:,:,:)
  real,dimension(numAttrib), intent(OUT) :: partAttribVec

  real, dimension(MDIM)   :: coord
  real                :: deltaCellInverseX, deltaCellInverseY, deltaCellInverseZ, xp, yp, zp
  real                :: hx, hy, hz, A, B, C, D, E, F, G, h1, h2
  integer             :: ip, jp, kp

  integer :: i, meshIndex
  real x,y,z
!-------------------------------------------------------------------------------


! Get a pointer to the grid data for the specified block, and obtain block
! size and coordinate information.

! Geometry issues
  coord(1) = bndBox(1,1)
  coord(2) = bndBox(1,2)
  coord(3) = bndBox(1,3)

  x=pos(IAXIS)
  y=pos(JAXIS)
  z=pos(KAXIS)

  if ( (pt_geometry /= CARTESIAN) .and. &
       (.not. ((pt_geometry == CYLINDRICAL) .and. (NDIM == 2))) ) &
    call Driver_abortFlash ("pt_mapFromMeshQuadratic:  unsupported geometry!")

!-------------------------------------------------------------------------------

! x-coord
  deltaCellInverseX = 1. / deltaCell(1)                ! inverse of deltaCell
  xp   = (x - coord(1)) * deltaCellInverseX     ! offset of particle from block edge
  ip   = floor(xp) + 1 + NGUARD    ! actual cell index (including guards)
  hx   = deltaCell(1) * (modulo(xp,1.) - 0.5)  ! remainder of something...??


#ifdef DEBUG_GRIDPARTICLES
  if(ip >= NXB+2*NGUARD .or. ip < 0) then   !!16 was NXB + 2NGUARD
     print '(4F8.2)', bndBox(1,1), bndBox(2,1), bndBox(1,2), bndBox(2,2)
     print *, "x, y = ", x, y
     !print *, "newBlockID = ", blockID
     print *, "xp = ", xp
     print *, "ip = ", ip
     call Driver_abortFlash("we have a problem in pt_mapFromMeshQuadratic")
  end if
#endif

! y-coord
  if (NDIM >= 2) then
    deltaCellInverseY = 1. / deltaCell(2)
    yp   = (y - coord(2)) * deltaCellInverseY
    jp   = floor(yp) + 1 + NGUARD
    hy   = deltaCell(2) * (modulo(yp,1.) - 0.5)
  else
    jp = 1
    hy = 0.
  endif


#ifdef DEBUG_GRIDPARTICLES
  if(jp >= NXB+2*NGUARD .or. jp < 0) then
     print '(4F8.2)', bndBox(1,1), bndBox(2,1), bndBox(1,2), bndBox(2,2)
     print *, "x, y = ", x, y
     !print *, "newBlockID = ", blockID
     print *, "yp = ", yp
     print *, "jp = ", jp
     call Driver_abortFlash("we have a problem in pt_mapFromMeshQuadratic")
  end if
#endif

! z-coord
  if (NDIM == 3) then
    deltaCellInverseZ = 1. / deltaCell(3)
    zp   = (z - coord(3)) * deltaCellInverseZ
    kp   = floor(zp) + 1 + NGUARD
    hz   = deltaCell(3) * (modulo(zp,1.) - 0.5)
  else
    kp = 1
    hz = 0.
  endif

!-------------------------------------------------------------------------------

  if (NDIM == 1) then
     
     do i = 1,numAttrib
        meshIndex=attrib(GRID_DS_IND,i)
        A = (    -solnVec(meshIndex,ip-1,jp  ,kp  ) &
             +26.*solnVec(meshIndex,ip  ,jp  ,kp  ) &
             -solnVec(meshIndex,ip+1,jp  ,kp  )) / 24.
        B = (    -solnVec(meshIndex,ip-1,jp  ,kp  ) &
             +solnVec(meshIndex,ip+1,jp  ,kp  )) * 0.5*deltaCellInverseX
        C = (     solnVec(meshIndex,ip-1,jp  ,kp  ) &
             -2.*solnVec(meshIndex,ip  ,jp  ,kp  ) &
             +solnVec(meshIndex,ip+1,jp  ,kp  )) * 0.5*deltaCellInverseX**2
        
        partAttribVec(i) = A + B*hx + C*hx**2
     end do
     
  endif
  
  !-------------------------------------------------------------------------------
  
  if (NDIM == 2) then
     
     if (pt_geometry == CARTESIAN) then
        do i = 1,numAttrib
           meshIndex = attrib(GRID_DS_IND,i)
           A = (    -solnVec(meshIndex,ip  ,jp-1,kp  ) &
                -solnVec(meshIndex,ip-1,jp  ,kp  ) &
                +28.*solnVec(meshIndex,ip  ,jp  ,kp  ) &
                -solnVec(meshIndex,ip+1,jp  ,kp  ) &
                -solnVec(meshIndex,ip  ,jp+1,kp  )) / 24.
           B = (    -solnVec(meshIndex,ip-1,jp  ,kp  ) &
                +solnVec(meshIndex,ip+1,jp  ,kp  )) * 0.5*deltaCellInverseX
           C = (     solnVec(meshIndex,ip-1,jp  ,kp  ) &
                -2.*solnVec(meshIndex,ip  ,jp  ,kp  ) &
                +solnVec(meshIndex,ip+1,jp  ,kp  )) * 0.5*deltaCellInverseX**2
           D = (    -solnVec(meshIndex,ip  ,jp-1,kp  ) &
                +solnVec(meshIndex,ip  ,jp+1,kp  )) * 0.5*deltaCellInverseY
           E = (     solnVec(meshIndex,ip  ,jp-1,kp  ) &
                -2.*solnVec(meshIndex,ip  ,jp  ,kp  ) &
                +solnVec(meshIndex,ip  ,jp+1,kp  )) * 0.5*deltaCellInverseY**2
           
           partAttribVec(i) = A + B*hx + C*hx**2 + D*hy + E*hy**2
        end do
     else   ! pt_geometry == CYLINDRICAL

        do i = 1,numAttrib

           meshIndex=attrib(GRID_DS_IND,i)
           h1 = ip - NGUARD - 0.5 + coord(1)*deltaCellInverseX         ! r_i / dr
           h2 = 4.*h1*h1 - 5.
           
           A = (                                  -solnVec(meshIndex,ip  ,jp-1,kp  ) &
                -(h1-1)/h1*solnVec(meshIndex,ip-1,jp  ,kp  ) &
                +28.*solnVec(meshIndex,ip  ,jp  ,kp  ) &
                -(h1+1)/h1*solnVec(meshIndex,ip+1,jp  ,kp  ) &
                -solnVec(meshIndex,ip  ,jp+1,kp  )) / 24.
           B = (               -(7+6*h1)*(h1-1)/h2*solnVec(meshIndex,ip-1,jp  ,kp  ) &
                +2*h1/h2*solnVec(meshIndex,ip  ,jp  ,kp  ) &
                +(-7+6*h1)*(1+h1)/h2*solnVec(meshIndex,ip+1,jp  ,kp  ))*deltaCellInverseX/3.
           C = ( (12*h1*h1+12*h1-1)*(h1-1)/(h1*h2)*solnVec(meshIndex,ip-1,jp  ,kp  ) &
                -2*(12*h1*h1-13)/h2*solnVec(meshIndex,ip  ,jp  ,kp  ) &
                +(12*h1*h1-12*h1-1)*(1+h1)/(h1*h2)*solnVec(meshIndex,ip+1,jp  ,kp  ))*deltaCellInverseX**2/6.
           D = (    -solnVec(meshIndex,ip  ,jp-1,kp  ) &
                +solnVec(meshIndex,ip  ,jp+1,kp  )) * 0.5*deltaCellInverseY
           E = (     solnVec(meshIndex,ip  ,jp-1,kp  ) &
                -2.*solnVec(meshIndex,ip  ,jp  ,kp  ) &
                +solnVec(meshIndex,ip  ,jp+1,kp  )) * 0.5*deltaCellInverseY**2
           
           partAttribVec(i) = A + B*hx + C*hx**2 + D*hy + E*hy**2
        end do
     endif

  endif

!-------------------------------------------------------------------------------

  if (NDIM == 3) then
     do i = 1,numAttrib
        meshIndex=attrib(GRID_DS_IND,i)
        A = (    -solnVec(meshIndex,ip  ,jp  ,kp-1) &
             -solnVec(meshIndex,ip  ,jp-1,kp  ) &
             -solnVec(meshIndex,ip-1,jp  ,kp  ) &
             +30.*solnVec(meshIndex,ip  ,jp  ,kp  ) &
             -solnVec(meshIndex,ip+1,jp  ,kp  ) &
             -solnVec(meshIndex,ip  ,jp+1,kp  ) &
             -solnVec(meshIndex,ip  ,jp  ,kp+1)) / 24.
        B = (    -solnVec(meshIndex,ip-1,jp  ,kp  ) &
             +solnVec(meshIndex,ip+1,jp  ,kp  )) * 0.5*deltaCellInverseX
        C = (     solnVec(meshIndex,ip-1,jp  ,kp  ) &
             -2.*solnVec(meshIndex,ip  ,jp  ,kp  ) &
             +solnVec(meshIndex,ip+1,jp  ,kp  )) * 0.5*deltaCellInverseX**2
        D = (    -solnVec(meshIndex,ip  ,jp-1,kp  ) &
             +solnVec(meshIndex,ip  ,jp+1,kp  )) * 0.5*deltaCellInverseY
        E = (     solnVec(meshIndex,ip  ,jp-1,kp  ) &
             -2.*solnVec(meshIndex,ip  ,jp  ,kp  ) &
             +solnVec(meshIndex,ip  ,jp+1,kp  )) * 0.5*deltaCellInverseY**2
        F = (    -solnVec(meshIndex,ip  ,jp  ,kp-1) &
             +solnVec(meshIndex,ip  ,jp  ,kp+1)) * 0.5*deltaCellInverseZ
        G = (     solnVec(meshIndex,ip  ,jp  ,kp-1) &
             -2.*solnVec(meshIndex,ip  ,jp  ,kp  ) &
             +solnVec(meshIndex,ip  ,jp  ,kp+1)) * 0.5*deltaCellInverseZ**2
        
        partAttribVec(i) = A + B*hx + C*hx**2 + D*hy + E*hy**2 + F*hz + G*hz**2
     end do
  endif
  
  !--------------------------------------------------------------------------
  
  return
  
end subroutine pt_mapFromMeshQuadratic

!===============================================================================

