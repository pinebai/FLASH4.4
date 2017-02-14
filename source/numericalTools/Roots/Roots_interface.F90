!!****h* source/numericalTools/Roots/Roots_interface
!!
!! NAME
!!
!!  Roots_interface
!!
!! SYNOPSIS
!!
!!  use Roots_interface
!!
!! DESCRIPTION
!!
!!  This is the header file for the Roots unit that defines its public interfaces.
!!
!!***

Module Roots_interface

  interface
    subroutine Roots_finalize ()
    end subroutine Roots_finalize
  end interface

  interface
    subroutine Roots_init ()
    end subroutine Roots_init
  end interface

  interface
    subroutine Roots_x2Polynomial (q1,q0,  nReal,root)
      real   , intent (in)  :: q1, q0
      integer, intent (out) :: nReal
      real   , intent (out) :: root (1:2,1:2)
    end subroutine Roots_x2Polynomial
  end interface

  interface
    subroutine Roots_x3Polynomial (c2,c1,c0,nReal,root,printInfo,printUnit)
      real   ,           intent (in)  :: c2, c1, c0
      integer,           intent (out) :: nReal
      real   ,           intent (out) :: root (1:3,1:2)
      logical, optional, intent (in)  :: printInfo
      integer, optional, intent (in)  :: printUnit
    end subroutine Roots_x3Polynomial
  end interface

  interface
    subroutine Roots_x4Polynomial (q3,q2,q1,q0,nReal,root,printInfo,printUnit)
      real   ,           intent (in)  :: q3, q2, q1, q0
      integer,           intent (out) :: nReal
      real   ,           intent (out) :: root (1:4,1:2)
      logical, optional, intent (in)  :: printInfo
      integer, optional, intent (in)  :: printUnit
    end subroutine Roots_x4Polynomial
  end interface

end Module Roots_interface
