!!****if* source/physics/Eos/EosMain/Eos_guardCells
!!
!! NAME
!!
!!  Eos_guardCells
!!
!! SYNOPSIS
!!
!!  call Eos_guardCells(integer(IN)  :: eosmode,
!!                      integer(IN)  :: blockid,
!!                      logical(IN)  :: corners,
!!             optional,integer(IN)  :: layers(MDIM),
!!             optional,logical(IN)  :: skipSrl)
!!
!! DESCRIPTION
!!
!!  Another layer of wrapping around Eos_wrapped calls, provided
!!  as a convenience to make it easy to apply the EOS to guard cells only.
!!  
!! ARGUMENTS
!!
!!  eosmode : determines which variables are used as Eos input variables.
!!            The valid values are MODE_DENS_EI (where density and internal
!!            energy are inputs), MODE_DENS_PRES (density and pressure as inputs)
!!            MODE_DENS_TEMP (density and temperature are inputs).
!!            These quantities are defined in constants.h.
!!            The argument is passed unchanged and unexamined to Eos_wrapped calls.
!!
!!  blockid : ID of block in current processor
!!
!!  corners : indicates whether Eos should be called on corner
!!            guard cells (i.e., diagonal guard cells)
!!  layers  : the number of guard cells to be included along each dimension
!!  skipSrl : whether to skip guard cell regions that are coming from
!!            neighboring blocks at the same refinement (or from boundary
!!            conditions) and thus have not undergone interpolation or
!!            restrictions.
!!
!! SEE ALSO
!!  Eos_wrapped
!!
!!***

subroutine Eos_guardCells(eosMode, blockID,corners,layers,skipSrl)

#include "Flash.h"
#include "constants.h"

  use Grid_interface, ONLY: Grid_getBlkIndexLimits, Grid_getBlkNeighLevels
  use Eos_interface,  ONLY : Eos_wrapped

  implicit none
  integer,intent(IN) :: eosMode,blockID
  logical,intent(IN) :: corners
  integer,dimension(MDIM),optional, intent(IN) :: layers
  logical,optional, intent(IN) :: skipSrl

  integer, dimension(2,MDIM) :: blkLimits,blkLimitsGC,eosRange
  integer,dimension(MDIM,2) :: nlayers
  integer :: neighLev(-1:1, -K2D:K2D , -K3D:K3D)
  integer :: myRefine
  logical :: skippingSrl

  if(present(skipSrl)) then
     skippingSrl = skipSrl
  else
     skippingSrl = .FALSE.
  end if

#ifdef FLASH_GRID_UG
  ! Nothing to do for a uniform Grid if skipping - we can take a shortcut, RETURN immediately.
  if (skippingSrl) then
     return
  end if
#endif

  call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)

  if(present(layers)) then
     nlayers(:,LOW)=layers
     nlayers(:,HIGH)=layers
  else
     nlayers(1:MDIM,LOW) =blkLimits(LOW,1:MDIM)   -blkLimitsGC(LOW,1:MDIM)
     nlayers(1:MDIM,HIGH)=blkLimitsGC(HIGH,1:MDIM)-blkLimits(HIGH,1:MDIM)
  end if

  if (skippingSrl) then
     ! Test whether neighbors are at same refinement level, and
     ! set the appropriate elements of nlayers to 0 when applicable.
     call Grid_getBlkNeighLevels(blockID,neighLev)
     myRefine = neighLev(0,0,0)
     if (.NOT.corners) then
        if (neighLev(-1,0,0) == myRefine) nlayers(IAXIS,LOW)  = 0
        if (neighLev( 1,0,0) == myRefine) nlayers(IAXIS,HIGH) = 0
#       if NDIM > 1
           if (neighLev(0,-1,0) == myRefine) nlayers(JAXIS,LOW)  = 0
           if (neighLev(0, 1,0) == myRefine) nlayers(JAXIS,HIGH) = 0
#       endif
#       if NDIM > 2
           if (neighLev(0,0,-1) == myRefine) nlayers(KAXIS,LOW)  = 0
           if (neighLev(0,0, 1) == myRefine) nlayers(KAXIS,HIGH) = 0
#       endif
     else                       !not corners
        if (ALL(neighLev(-1,:,:) == myRefine)) nlayers(IAXIS,LOW)  = 0
        if (ALL(neighLev( 1,:,:) == myRefine)) nlayers(IAXIS,HIGH) = 0
#       if NDIM > 1
           if (ALL(neighLev(:,-1,:) == myRefine)) nlayers(JAXIS,LOW)  = 0
           if (ALL(neighLev(:, 1,:) == myRefine)) nlayers(JAXIS,HIGH) = 0
#       endif
#       if NDIM > 2
           if (ALL(neighLev(:,:,-1) == myRefine)) nlayers(KAXIS,LOW)  = 0
           if (ALL(neighLev(:,:, 1) == myRefine)) nlayers(KAXIS,HIGH) = 0
#       endif
     end if
  end if

  eosRange = blkLimits

  if (NDIM>1 .AND. corners .AND. skippingSrl) then
     call complexSkipping()
     return
  end if


  eosRange(LOW,IAXIS) = blkLimits(LOW,IAXIS)-nlayers(IAXIS,LOW)
  eosRange(HIGH,IAXIS) = blkLimits(LOW,IAXIS)-1
  call Eos_wrapped(eosMode,eosRange,blockID)
  eosRange(LOW,IAXIS) = blkLimits(HIGH,IAXIS)+1
  eosRange(HIGH,IAXIS) = blkLimits(HIGH,IAXIS)+nlayers(IAXIS,HIGH)
  call Eos_wrapped(eosMode,eosRange,blockID)

# if NDIM > 1
     if (corners) then
        eosRange(LOW,IAXIS) = blkLimits(LOW,IAXIS)-nlayers(IAXIS,LOW)
        eosRange(HIGH,IAXIS) = blkLimits(HIGH,IAXIS)+nlayers(IAXIS,HIGH)
     else
        eosRange(:,IAXIS) = blkLimits(:,IAXIS)
     end if
     eosRange(LOW,JAXIS) = blkLimits(LOW,JAXIS)-nlayers(JAXIS,LOW)
     eosRange(HIGH,JAXIS) = blkLimits(LOW,JAXIS)-1
     call Eos_wrapped(eosMode,eosRange,blockID)
     eosRange(LOW,JAXIS) = blkLimits(HIGH,JAXIS)+1
     eosRange(HIGH,JAXIS) = blkLimits(HIGH,JAXIS)+nlayers(JAXIS,HIGH)
     call Eos_wrapped(eosMode,eosRange,blockID)
# endif

# if NDIM > 2
     if (corners) then
        eosRange(LOW,JAXIS) = blkLimits(LOW,JAXIS)-nlayers(JAXIS,LOW)
        eosRange(HIGH,JAXIS) = blkLimits(HIGH,JAXIS)+nlayers(JAXIS,HIGH)
     else
        eosRange(:,JAXIS) = blkLimits(:,JAXIS)
     end if
     eosRange(LOW,KAXIS) = blkLimits(LOW,KAXIS)-nlayers(KAXIS,LOW)
     eosRange(HIGH,KAXIS) = blkLimits(LOW,KAXIS)-1
     call Eos_wrapped(eosMode,eosRange,blockID)
     eosRange(LOW,KAXIS) = blkLimits(HIGH,KAXIS)+1
     eosRange(HIGH,KAXIS) = blkLimits(HIGH,KAXIS)+nlayers(KAXIS,HIGH)
     call Eos_wrapped(eosMode,eosRange,blockID)
# endif

contains
  subroutine complexSkipping()
    use Eos_data,       ONLY : eos_meshMe

    logical :: done(-1:1, -K2D:K2D , -K3D:K3D)
    integer :: i,j,k, il,jl,kl, ir,jr,kr
    integer :: d,s

    integer,dimension(2,-1:1,MDIM) :: surr2lim
    integer :: numEosCalls

    do s=-1,1
       surr2lim(:,s,:) = blkLimits(:,:)
    end do
    do d=1,NDIM
       surr2lim(LOW ,-1,d) = blkLimits(LOW,d)  - nlayers(d,LOW)
       surr2lim(HIGH,-1,d) = blkLimits(LOW,d)  - 1
       surr2lim(LOW , 1,d) = blkLimits(HIGH,d) + 1
       surr2lim(HIGH, 1,d) = blkLimits(HIGH,d) + nlayers(d,HIGH)
    end do

    done(:,:,:) = .FALSE.
    done(0,0,0) = .TRUE.
    where (neighLev(:,:,:) == myRefine) done(:,:,:) = .TRUE.
    if (nlayers(IAXIS,LOW)  == 0) done(-1,:,:) = .TRUE.
    if (nlayers(IAXIS,HIGH) == 0) done( 1,:,:) = .TRUE.
#   if NDIM > 1
       if (nlayers(JAXIS,LOW)  == 0) done(:,-1,:) = .TRUE.
       if (nlayers(JAXIS,HIGH) == 0) done(:, 1,:) = .TRUE.
#   endif
#   if NDIM > 2
       if (nlayers(KAXIS,LOW)  == 0) done(:,:,-1) = .TRUE.
       if (nlayers(KAXIS,HIGH) == 0) done(:,:, 1) = .TRUE.
#   endif

    numEosCalls = 0

#   if NDIM > 2
       do k=-1,1,2
          if (.NOT.done(0,0,k)) then
             eosRange(:,IAXIS:JAXIS) = blkLimits(:,IAXIS:JAXIS)
             jl = 0; jr = 0
             if (.NOT.done(0,-1,k)) then
                eosRange(LOW,JAXIS) = blkLimits(LOW,JAXIS)-nlayers(JAXIS,LOW)
                jl = -1
             end if
             if (.NOT.done(0, 1,k)) then
                eosRange(HIGH,JAXIS) = blkLimits(HIGH,JAXIS)+nlayers(JAXIS,HIGH)
                jr = 1
             end if
             il = 0; ir = 0
             if (ALL(.NOT.done(-1,jl:jr,k))) then
                eosRange(LOW,IAXIS) = blkLimits(LOW,IAXIS)-nlayers(IAXIS,LOW)
                il = -1
             end if
             if (ALL(.NOT.done( 1,jl:jr,k))) then
                eosRange(HIGH,IAXIS) = blkLimits(HIGH,IAXIS)+nlayers(IAXIS,HIGH)
                ir =  1
             end if
             eosRange(:,KAXIS) = surr2lim(:,k,KAXIS)
             call Eos_wrapped(eosMode,eosRange,blockID); numEosCalls = numEosCalls+1
             done(il:ir,jl:jr,k) = .TRUE.
          end if
       end do
#   endif

    do j=-1,1,2
       if (.NOT.done(0,j,0)) then
          eosRange(:,IAXIS) = blkLimits(:,IAXIS)
          kl = 0; kr = 0
#         if NDIM > 2
             eosRange(:,KAXIS) = blkLimits(:,KAXIS)
             if (.NOT.done(0,j,-1)) then
                eosRange(LOW,KAXIS) = blkLimits(LOW,KAXIS)-nlayers(KAXIS,LOW)
                kl = -1
             end if
             if (.NOT.done(0,j, 1)) then
                eosRange(HIGH,KAXIS) = blkLimits(HIGH,KAXIS)+nlayers(KAXIS,HIGH)
                kr = 1
             end if
#         endif
          il = 0; ir = 0
          if (ALL(.NOT.done(-1,j,kl:kr))) then
             eosRange(LOW,IAXIS) = blkLimits(LOW,IAXIS)-nlayers(IAXIS,LOW)
             il = -1
          end if
          if (ALL(.NOT.done( 1,j,kl:kr))) then
             eosRange(HIGH,IAXIS) = blkLimits(HIGH,IAXIS)+nlayers(IAXIS,HIGH)
             ir =  1
          end if
          eosRange(:,JAXIS) = surr2lim(:,j,JAXIS)
          call Eos_wrapped(eosMode,eosRange,blockID); numEosCalls = numEosCalls+1
          done(il:ir,j,kl:kr) = .TRUE.
       end if
    end do

    do i=-1,1,2
       if (.NOT.done(i,0,0)) then
          eosRange(:,JAXIS) = blkLimits(:,JAXIS)
          kl = 0; kr = 0
#         if NDIM > 2
             eosRange(:,KAXIS) = blkLimits(:,KAXIS)
             if (.NOT.done(i,0,-1)) then
                eosRange(LOW,KAXIS) = blkLimits(LOW,KAXIS)-nlayers(KAXIS,LOW)
                kl = -1
             end if
             if (.NOT.done(i,0, 1)) then
                eosRange(HIGH,KAXIS) = blkLimits(HIGH,KAXIS)+nlayers(KAXIS,HIGH)
                kr = 1
             end if
#         endif
          jl = 0; jr = 0
#         if NDIM > 1
          if (ALL(.NOT.done(i,-1,kl:kr))) then
             eosRange(LOW,JAXIS) = blkLimits(LOW,JAXIS)-nlayers(JAXIS,LOW)
             jl = -1
          end if
          if (ALL(.NOT.done(i, 1,kl:kr))) then
             eosRange(HIGH,JAXIS) = blkLimits(HIGH,JAXIS)+nlayers(JAXIS,HIGH)
             jr =  1
          end if
#         endif
          eosRange(:,IAXIS) = surr2lim(:,i,IAXIS)
          call Eos_wrapped(eosMode,eosRange,blockID); numEosCalls = numEosCalls+1
          done(i,jl:jr,kl:kr) = .TRUE.
       end if
    end do


#   if NDIM > 2
       do j=-1,1,2
          do i=-1,1,2
             if (.NOT.done(i,j,0)) then
                kl = 0; kr = 0
                eosRange(:,KAXIS) = blkLimits(:,KAXIS)
                if (.NOT.done(i,j,-1)) then
                   eosRange(LOW,KAXIS) = blkLimits(LOW,KAXIS)-nlayers(KAXIS,LOW)
                   kl = -1
                end if
                if (.NOT.done(i,j, 1)) then
                   eosRange(HIGH,KAXIS) = blkLimits(HIGH,KAXIS)+nlayers(KAXIS,HIGH)
                   kr = 1
                end if
                eosRange(:,IAXIS) = surr2lim(:,i,IAXIS)
                eosRange(:,JAXIS) = surr2lim(:,j,JAXIS)
                call Eos_wrapped(eosMode,eosRange,blockID); numEosCalls = numEosCalls+1
                done(i,j,kl:kr) = .TRUE.
             end if
          end do
       end do
#   endif



    if (ANY(.NOT.done(:,:,:))) then
       do k=-K3D,K3D
          do j=-1,1
             do i=-1,1
                if (.NOT.done(i,j,k)) then
                   eosRange(:,IAXIS) = surr2lim(:,i,IAXIS)
                   eosRange(:,JAXIS) = surr2lim(:,j,JAXIS)
                   eosRange(:,KAXIS) = surr2lim(:,j,KAXIS)
                   call Eos_wrapped(eosMode,eosRange,blockID); numEosCalls = numEosCalls+1
                end if
             end do
          end do
       end do
    end if

#ifdef DEBUG_EOS
!    if (numEosCalls>2*NDIM) then
       print*,'Eos_guardCells:',numEosCalls,' Eos_wr calls on block',blockID,' @',eos_meshMe
!       print*,neighLev
!    end if
#endif

  end subroutine complexSkipping

end subroutine Eos_guardCells
