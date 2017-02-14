module bittree
  implicit none
  
  interface
    subroutine amr_build_bittree
    end subroutine amr_build_bittree
  end interface
  
  interface
    subroutine amr_identify_block(procs, lev, ijk, proc, locblk)
      integer, intent(in) :: procs
      integer, intent(inout) :: lev
      integer, intent(inout) :: ijk(*)
      integer, intent(out) :: proc
      integer, intent(out) :: locblk
    end subroutine
  end interface
  
  interface
    subroutine amr_verify_bittree
    end subroutine
  end interface
  
  interface
    function bittree_initialized() result(yep) bind(c,name='bittree_initialized')
      use iso_c_binding, only: c_bool
      logical(c_bool) :: yep
    end function
  end interface
  
  interface
    subroutine bittree_init(ndim, top, topmask) bind(c,name='bittree_init')
      use iso_c_binding, only: c_bool, c_int
      integer(c_int), intent(in) :: ndim
      integer(c_int), intent(in) :: top(*)
      logical(c_bool), intent(in) :: topmask(*)
    end subroutine
  end interface
  
  interface
    subroutine bittree_block_count(count) bind(c,name='bittree_block_count')
      use iso_c_binding, only: c_int
      integer(c_int), intent(out) :: count
    end subroutine
  end interface
  
  interface
    subroutine bittree_morton(lev, ijk, mort) bind(c,name='bittree_morton')
      use iso_c_binding, only: c_int
      integer(c_int), intent(inout) :: lev
      integer(c_int), intent(inout) :: ijk(*)
      integer(c_int), intent(out) :: mort
    end subroutine
  end interface
  
  interface
    subroutine bittree_refine_init() bind(c,name='bittree_refine_init')
    end subroutine
  end interface
  
  interface
    subroutine bittree_refine_mark(lev, ijk) bind(c,name='bittree_refine_mark')
      use iso_c_binding, only: c_int
      integer(c_int), intent(in) :: lev
      integer(c_int), intent(in) :: ijk(*)
    end subroutine
  end interface
  
  interface
    subroutine bittree_refine_apply(comm) bind(c,name='bittree_refine_apply')
      use iso_c_binding, only: c_int
      integer(c_int), intent(in) :: comm
    end subroutine
  end interface
end module
