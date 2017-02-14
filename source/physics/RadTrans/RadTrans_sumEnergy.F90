!!****f* source/physics/RadTrans/RadTrans_sumEnergy
!!
!!  NAME 
!!
!!  RadTrans_sumEnergy
!!
!!  SYNOPSIS
!!
!!  call RadTrans_sumEnergy( integer(IN) :: ivar,
!!                           integer(IN) :: nblk,
!!                           integer(IN) :: blklst(nblk),
!!
!!  DESCRIPTION 
!!
!!    This subroutine is useful when mesh replication is active
!!    (meshCopyCount > 1). It takes an unk variable and adds it across
!!    all the meshes. For example, if each mesh had computed ERAD
!!    separately using its own energy groups, this subroutine could be
!!    used to add all of the ERADs to compute the total radiation
!!    energy.
!!
!! ARGUMENTS
!!
!!   ivar   : the unk variable to add
!!   blklst : the list of blocks to cover
!!   nblk   : the number of blocks
!!
!!***
subroutine RadTrans_sumEnergy(ivar, nblk, blklst)
  implicit none
  integer, intent(in) :: ivar
  integer, intent(in) :: nblk
  integer, intent(in) :: blklst(nblk)
  ! Stub Implementation
end subroutine RadTrans_sumEnergy
     
