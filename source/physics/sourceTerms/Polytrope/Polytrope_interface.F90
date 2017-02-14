!!****h* source/physics/sourceTerms/Polytrope/Polytrope_interface
!!
!! This is the header file for the Polytrope module that defines its
!! public interfaces.
!!
!! WRITTEN BY
!!   Christoph Federrath 2007
!!
!!***

Module Polytrope_interface

  interface Polytrope
    subroutine Polytrope(blockCount,blockList,dt)
      integer, intent(IN) :: blockCount
      integer, dimension(blockCount), intent(IN) :: blockList
      real, intent(IN) :: dt
    end subroutine Polytrope
  end interface

  interface Polytrope_finalize
    subroutine Polytrope_finalize()
    end subroutine Polytrope_finalize
  end interface

  interface Polytrope_init
    subroutine Polytrope_init()
    end subroutine Polytrope_init
  end interface

end Module Polytrope_interface


