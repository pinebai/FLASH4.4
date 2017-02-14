!!****if* source/physics/sourceTerms/Flame/FlameSpeed/BuoyancyCompensation/fl_fsAtwood
!!
!! NAME
!!
!!  fl_fsAtwood
!!
!! SYNOPSIS
!!
!!  call fl_fsAtwood(real, dimension(:,:,:,:), pointer  :: solndata,
!!                   real, dimension(:,:,:)(out) :: atwood,
!!                   real, dimension(:,:,:)(in) :: dens_u,
!!                   integer, dimension(LOW:HIGH,MDIM)(in) :: cl)
!!
!! DESCRIPTION
!!
!!
!! ARGUMENTS
!!
!!   solndata : 
!!
!!   atwood : 
!!
!!   dens_u : 
!!
!!   cl : 
!!
!! AUTOGENROBODOC
!!
!!
!!***


! Copied from trunk Flame code - KW 2013

subroutine fl_fsAtwood(solnData, atwood, dens_u, cl)

#include "constants.h"
#include "Flash.h"
  implicit none

  real, dimension(:,:,:,:), pointer :: solnData
  real, dimension(:,:,:), intent(out) :: atwood
  real, dimension(:,:,:), intent(in)  :: dens_u
  integer, dimension(LOW:HIGH,MDIM), intent(in)    :: cl

  logical,parameter :: fl_approxAtwood = .FALSE.

#ifdef FIXEDBLOCKSIZE 
!!$  real, dimension(GRID_ILO_GC:GRID_IHI_GC,&
!!$                  GRID_JLO_GC:GRID_JHI_GC,&
!!$                  GRID_KLO_GC:GRID_KHI_GC) :: temp
  real, dimension(cl(LOW,1):cl(HIGH,1),cl(LOW,2):cl(HIGH,2),cl(LOW,3):cl(HIGH,3)) :: &
       factor,phi,qburn,dens,temp
#else
  real, allocatable, dimension(:,:,:) :: factor,phi,qburn,dens,temp
#endif

#ifndef FIXEDBLOCKSIZE
     allocate(factor(cl(LOW,1):cl(HIGH,1),cl(LOW,2):cl(HIGH,2),cl(LOW,3):cl(HIGH,3)))
     allocate(phi(cl(LOW,1):cl(HIGH,1),cl(LOW,2):cl(HIGH,2),cl(LOW,3):cl(HIGH,3)))
     allocate(qburn(cl(LOW,1):cl(HIGH,1),cl(LOW,2):cl(HIGH,2),cl(LOW,3):cl(HIGH,3)))
     allocate(dens(cl(LOW,1):cl(HIGH,1),cl(LOW,2):cl(HIGH,2),cl(LOW,3):cl(HIGH,3)))
     allocate(temp(cl(LOW,1):cl(HIGH,1),cl(LOW,2):cl(HIGH,2),cl(LOW,3):cl(HIGH,3)))
#endif

  dens = solnData(DENS_VAR,cl(LOW,1):cl(HIGH,1),cl(LOW,2):cl(HIGH,2),cl(LOW,3):cl(HIGH,3))
  temp = solnData(TEMP_VAR,cl(LOW,1):cl(HIGH,1),cl(LOW,2):cl(HIGH,2),cl(LOW,3):cl(HIGH,3))
#ifdef RPV1_MSCALAR
  phi  = solnData(RPV1_MSCALAR, cl(LOW,1):cl(HIGH,1),cl(LOW,2):cl(HIGH,2),cl(LOW,3):cl(HIGH,3))
#else
  phi  = solnData(FLAM_MSCALAR, cl(LOW,1):cl(HIGH,1),cl(LOW,2):cl(HIGH,2),cl(LOW,3):cl(HIGH,3))
#endif

!.. laminar

  call heatReleaseBlock(qburn)
 
  factor = solnData(GAME_VAR,cl(LOW,1):cl(HIGH,1),cl(LOW,2):cl(HIGH,2),cl(LOW,3):cl(HIGH,3))
  factor = (factor - 1.e0)/factor
  factor = factor*qburn*dens/solnData(PRES_VAR,cl(LOW,1):cl(HIGH,1),cl(LOW,2):cl(HIGH,2),cl(LOW,3):cl(HIGH,3))

  if (fl_approxAtwood) then
     atwood(cl(LOW,1):cl(HIGH,1),cl(LOW,2):cl(HIGH,2),cl(LOW,3):cl(HIGH,3)) = 0.5e0*factor
  else
     factor = 1.e0 + max(0.e0, factor/(1.e0 - factor*phi))
     !.atwood number
     atwood(cl(LOW,1):cl(HIGH,1),cl(LOW,2):cl(HIGH,2),cl(LOW,3):cl(HIGH,3)) = (factor - 1.e0)/(factor + 1.e0)
  endif

#ifndef FIXEDBLOCKSIZE
  deallocate(factor)
  deallocate(phi)
  deallocate(qburn)
  deallocate(dens)
  deallocate(temp)
#endif

contains
  subroutine heatReleaseBlock(q)

    use bn_paraData, ONLY : pbQc2mg, pbC12Init

    implicit none

    real, dimension(cl(LOW,1):cl(HIGH,1),cl(LOW,2):cl(HIGH,2),cl(LOW,3):cl(HIGH,3)),   intent(out) :: q

    q = pbQc2mg * pbC12Init

  end subroutine heatReleaseBlock

end subroutine
