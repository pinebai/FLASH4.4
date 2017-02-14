!!****if* source/physics/Gravity/GravityMain/Poisson/BHTree/Gravity_bhGetNodeStruct
!!
!! NAME
!!
!!  Gravity_bhGetNodeStruct
!!
!!
!! SYNOPSIS
!!
!!   call Gravity_bhGetNodeStruct(
!!                  integer(in)    :: im,
!!                  integer(in)    :: ix,
!!                  integer(in)    :: iy,
!!                  integer(in)    :: iz,
!!                  integer(inout) :: nsize,
!!                  integer(inout) :: bnsize
!!        )
!!
!! DESCRIPTION
!!
!!   Calculates structure (size, indeces of variables) of the 
!!   gravity section of the tree node.
!!
!! ARGUMENTS
!!
!!  im          : index of the mass in the node array 
!!  ix          : index x coord of the mass centre in the node array 
!!  iy          : index y coord of the mass centre in the node array 
!!  iz          : index z coord of the mass centre in the node array 
!!  nsize       : size of the node array (initially, before the gravity part is
!!                added; at return, after the gravity part is added)
!!  bnsize      : size of the bottom-most node array (initially, before 
!!                the gravity part is added; at return, after the 
!!                gravity part is added)
!!
!!
!!***

subroutine Gravity_bhGetNodeStruct(im, ix, iy, iz, nsize, bnsize)
  use Gravity_data, ONLY : grv_bhMAC, grv_bhUseRelAccErr, grv_bhNSIZE, &
    grv_bhBNSIZE, grv_bhNODE5, grv_bhIB2, grv_bhIB3, grv_bhIDMIN, grv_bhN5_B2, &
    grv_bhN5_DMIN, grv_bhN5_NONE, grv_bhIM, grv_bhIBM, grv_bhIX, grv_bhIY, grv_bhIZ, &
    useGravity, grv_bhPhysMACTW
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
  implicit none
  integer, intent(IN) :: im, ix, iy, iz
  integer, intent(INOUT) :: nsize, bnsize
  
  ! check whether gravity module is used
  call RuntimeParameters_get("useGravity", useGravity)
  if (.not. useGravity) return

  ! get some runtime params., this subroutine is called before Gravity_init
  call RuntimeParameters_get("gr_bhPhysMACTW", grv_bhPhysMACTW)
  call RuntimeParameters_get("grv_bhMAC", grv_bhMAC)
  call RuntimeParameters_get("grv_bhUseRelAccErr", grv_bhUseRelAccErr)

  ! set indeces of grid section quantities
  grv_bhIM  = im
  grv_bhIBM = im
  grv_bhIX  = ix
  grv_bhIY  = iy
  grv_bhIZ  = iz

  if (grv_bhPhysMACTW .and. .not. grv_bhUseRelAccErr) then
    ! set size of the gravity section
    grv_bhNSIZE  = 1
    grv_bhBNSIZE = 0

    ! set indeces of the additional quantities in the node
    if (grv_bhMAC == "SumSquare") then
      ! with SumSquare MAC, node(5) includes just B2 moment
      grv_bhNODE5 = grv_bhN5_B2
      grv_bhIB2   = nsize + 1
      grv_bhIDMIN = -1
      !grv_bhIB3   = GR_TREE_BASE_NSIZE + 2
    else
      ! with other MACs, node(5) includes minimum distance for the node acceptance
      ! grv_bhIB2 is also set, because DMIN is calculated from B2 and the code 
      ! for B2 calculation is shared
      grv_bhNODE5 = grv_bhN5_DMIN
      grv_bhIDMIN = nsize + 1
      grv_bhIB2   = grv_bhIDMIN
      !grv_bhIB3   = GR_TREE_BASE_NSIZE + 2
    endif
  else
    ! with relative acceleration error, nothing is added by gravity to the node
    grv_bhNSIZE  = 0
    grv_bhBNSIZE = 0
    grv_bhNODE5 = grv_bhN5_NONE
  endif

  ! total size of the node returned
  nsize  = nsize  + grv_bhNSIZE
  bnsize = bnsize + grv_bhBNSIZE
  grv_bhNSIZE  = nsize
  grv_bhBNSIZE = bnsize

  return
end subroutine Gravity_bhGetNodeStruct
