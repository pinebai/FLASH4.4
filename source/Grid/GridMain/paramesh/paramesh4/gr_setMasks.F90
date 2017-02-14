!!****if* source/Grid/GridMain/paramesh/paramesh4/gr_setMasks
!!
!! NAME
!!  gr_setMasks
!!
!! SYNOPSIS
!!
!!  gr_setMasks(integer(IN) :: gridDataStruct, 
!!                      integer(IN) :: maskSize,
!!                      logical(IN) :: mask(maskSize))
!!  
!! DESCRIPTION 
!!  
!!  This routine takes the mask supplied by the calling routine
!!  and sets the paramesh3 masks for cell centered and face centered
!!  variables as specified
!!
!!
!! ARGUMENTS 
!!  
!!
!!  gridDataStruct - integer constant, defined in "constants.h", 
!!                   indicating which grid data structure 
!!                   variable's guardcells to fill.
!!                   Paramesh has 5 data structures for grid variables, the first
!!                   four include all physical variables defined on the mesh. The 
!!                   fifth one includes a single variable.
!!
!!                   unk                cell centered, 
!!                   facex,facey,facez  face centered along i,j,k 
!!                                      direction respectively
!!                   work               cell centered, single variable.
!!                   
!!                   valid values of gridDataStruct are  
!!                   CENTER             unk only
!!                   WORK               work 
!!                   FACEX              facex
!!                   FACEY              facey
!!                   FACEZ              facez
!!                   FACES              facex,facey and facez
!!                   CENTER_FACES     unk,facex,facey,facez
!!
!!  maskSize - the size of the mask array. 
!! 
!!  mask -  It is a one-dimensional logical array 
!!          with indices correspoding to variables in the grid data
!!          structures. If a variable should have its guardcells filled,
!!          the corresponding element in "mask" is true, otherwise it is
!!          false.
!!          The mask is always ignored if the runtime parameter
!!          enableMaskedGCFill is set .FALSE.
!! 
!!
!! 
!!***


subroutine gr_setMasks(gridDataStruct,maskSize,mask)

#include "Flash.h"

  use Grid_data, ONLY : gr_enableMaskedGCFill

  !! This section disabled for now, may become useful with non-permanent gc support
!!$#ifdef FL_NON_PERMANENT_GUARDCELLS
!!$  use Grid_data, ONLY : gr_ccMask
!!$#if NFACE_VARS>0
!!$  use Grid_data, ONLY : gr_fcMask
!!$#endif
!!$#endif

  use Driver_interface, ONLY : Driver_abortFlash
  use physicaldata, ONLY : gcell_on_cc,gcell_on_fc, no_permanent_guardcells

  implicit none
#include "constants.h"
#include "Flash_mpi.h"

  integer, intent(in) :: gridDataStruct
  integer, intent(in) :: maskSize
  logical,dimension(maskSize),intent(in) :: mask
  
  integer :: nv,i

#ifdef DEBUG_GRID
  logical:: validDataStructure,validMaskSize

  validDataStructure = (gridDataStruct==CENTER).or.&
                       (gridDataStruct==FACES).or.&
                       (gridDataStruct==FACEX).or.&
                       (gridDataStruct==FACEY).or.&
                       (gridDataStruct==FACEZ).or.&
                       (gridDataStruct==WORK).or.&
                       (gridDataStruct==CENTER_FACES)
  if(.not.validDataStructure)then
     call Driver_abortFlash("GCfill: invalid data structure")
  end if
  select case(gridDataStruct) 
  case(CENTER)
     validMaskSize=(maskSize==NUNK_VARS)
  case(WORK)
     validMaskSize=.TRUE.    !but whatever the mask is, it will be ignored. 
  case(FACES)
     validMaskSize= ((maskSize.GE.NDIM*NFACE_VARS) .AND. (maskSize.LE.MDIM*NFACE_VARS))
  case(FACEX)
     validMaskSize=(maskSize==NFACE_VARS)
  case(FACEY)
     validMaskSize=(maskSize==NFACE_VARS)
  case(FACEZ)
     validMaskSize=(maskSize==NFACE_VARS)
  case(CENTER_FACES)
     validMaskSize= ((maskSize.GE.(NUNK_VARS+NDIM*NFACE_VARS)) .AND. (maskSize.LE.(NUNK_VARS+MDIM*NFACE_VARS)))
  end select
  if(.not.validMaskSize)then
     call Driver_abortFlash("GCfill : mask size mismatch")
  end if
#endif
  
  
  if (.NOT. gr_enableMaskedGCFill) return
  
  if((gridDataStruct==CENTER_FACES).or.(gridDataStruct==CENTER)) then
     gcell_on_cc(1:NUNK_VARS)=mask(1:NUNK_VARS)
  end if
  
#if NFACE_VARS>0 
  select case(gridDataStruct)
  case(FACEX)
     gcell_on_fc(IAXIS,1:NFACE_VARS)=mask(1:NFACE_VARS)
  case(FACEY)
     gcell_on_fc(JAXIS,1:NFACE_VARS)=mask(1:NFACE_VARS)
  case(FACEZ)
     gcell_on_fc(KAXIS,1:NFACE_VARS)=mask(1:NFACE_VARS)
  case(FACES)
     do i=1,NDIM
        gcell_on_fc(i,1:NFACE_VARS)=mask((i-1)*NFACE_VARS+1:i*NFACE_VARS)
     end do
     if (maskSize > NDIM*NFACE_VARS) then
        do i=NDIM+1,(maskSize)/(NFACE_VARS)
           gcell_on_fc(i,1:NFACE_VARS)=mask((i-1)*NFACE_VARS+1:i*NFACE_VARS)
        end do
        do i=maskSize/(NFACE_VARS)+1,(maskSize-1)/(NFACE_VARS)+1
           nv = mod(maskSize,NFACE_VARS)
           gcell_on_fc(i,1:nv)=mask((i-1)*NFACE_VARS+1:(i-1)*NFACE_VARS+nv)
        end do
     end if
  case(CENTER_FACES)
     do i=1,NDIM
        gcell_on_fc(i,1:NFACE_VARS)=mask(NUNK_VARS+(i-1)*NFACE_VARS+1:NUNK_VARS+i*NFACE_VARS)
     end do
     if (maskSize > NUNK_VARS+NDIM*NFACE_VARS) then
        do i=NDIM+1,(maskSize-NUNK_VARS)/(NFACE_VARS)
           gcell_on_fc(i,1:NFACE_VARS)=mask(NUNK_VARS+(i-1)*NFACE_VARS+1:NUNK_VARS+i*NFACE_VARS)
        end do
        do i=(maskSize-NUNK_VARS)/(NFACE_VARS)+1,(maskSize-NUNK_VARS-1)/(NFACE_VARS)+1
           nv = mod(maskSize-NUNK_VARS,NFACE_VARS)
           gcell_on_fc(i,1:nv)=mask(NUNK_VARS+(i-1)*NFACE_VARS+1:NUNK_VARS+(i-1)*NFACE_VARS+nv)
        end do
     end if
  end select
#endif

!! NOTE -- This section is currently redundant. It should be resurrected if we support
!!         non permanent guardcells.
!!$  if (no_permanent_guardcells) then
!!$     
!!$#ifdef FL_NON_PERMANENT_GUARDCELLS
!!$     gr_ccMask=gcell_on_cc
!!$#if NFACE_VARS >0
!!$     gr_fcMask=gcell_on_fc
!!$#endif
!!$#endif

end subroutine gr_setMasks
