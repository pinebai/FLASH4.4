!!****if* source/Grid/GridMain/UG/physicaldata
!!
!! NAME
!!  physicaldata
!!
!! SYNOPSIS
!!
!!  use physicaldata
!!
!! DESCRIPTION
!!  
!!  This is the data modules that defines data structures for the 
!!  physical variables on the discretized mesh for Uniform Grid. 
!!  The cell centered data structure is "unk", the three face centered
!!  ones are "facevarx", "facevary" and "facevarz". There is also
!!  a data strucure to provide global domain scope scratch space called
!!  "scratch"
!!
!!***

! solnData depends on the ordering on unk
!!REORDER(5): unk, facevar[xyz]


module physicalData
#include "Flash.h"

#ifdef FIXEDBLOCKSIZE
  real, save, target :: &
       unk(UNK_VARS_BEGIN:UNK_VARS_END,&
       GRID_ILO_GC:GRID_IHI_GC, GRID_JLO_GC:GRID_JHI_GC,&
       GRID_KLO_GC:GRID_KHI_GC,1) 

#if(NFACE_VARS>0)

  real, save, target :: &
       facevarx( NFACE_VARS,&
       GRID_ILO_GC:GRID_IHI_GC+1, GRID_JLO_GC:GRID_JHI_GC,&
       GRID_KLO_GC:GRID_KHI_GC,1) 
  real, save, target :: &
       facevary( NFACE_VARS,&
       GRID_ILO_GC:GRID_IHI_GC, GRID_JLO_GC:GRID_JHI_GC+K2D,&
       GRID_KLO_GC:GRID_KHI_GC,1) 
  real, save, target :: &
       facevarz( NFACE_VARS,&
       GRID_ILO_GC:GRID_IHI_GC, GRID_JLO_GC:GRID_JHI_GC,&
       GRID_KLO_GC:GRID_KHI_GC+K3D,1) 

#else
  real, save,target,dimension(1,1,1,1,1)::facevarx,facevary,facevarz
#endif

#else
  real, save, target, allocatable :: unk(:,:,:,:,:)
  real, save, target, allocatable :: facevarx(:,:,:,:,:)
  real, save, target, allocatable :: facevary(:,:,:,:,:)
  real, save, target, allocatable :: facevarz(:,:,:,:,:)
#endif

end module physicalData
