!!****h* source/numericalTools/Interpolate/Interpolate_interface
!!
!! NAME
!!
!!  Interpolate_interface
!!
!! SYNOPSIS
!!
!!   use Interpolate_interface
!!
!! DESCRIPTION
!!
!!  This is the module needed for defining the interpolation unit interface.
!!
!!***

Module Interpolate_interface

  interface
    subroutine Interpolate_cubic1Dcoeffs (numberOfLines, a)
      integer, intent (in)    :: numberOfLines
      real,    intent (inout) :: a (1:4,*)
    end subroutine Interpolate_cubic1Dcoeffs
  end interface

  interface
    function Interpolate_cubic1DFd1d2 (a,x)
      real, intent (in) :: a (1:4)
      real, intent (in) :: x
      real :: Interpolate_cubic1DFd1d2 (1:3)
    end function Interpolate_cubic1DFd1d2
  end interface

  interface
    function Interpolate_cubic1DFd1 (a,x)
      real, intent (in) :: a (1:4)
      real, intent (in) :: x
      real :: Interpolate_cubic1DFd1 (1:2)
    end function Interpolate_cubic1DFd1
  end interface

  interface
    real function Interpolate_cubic1DF (a,x)
      real, intent (in) :: a (1:4)
      real, intent (in) :: x
    end function Interpolate_cubic1DF
  end interface

  interface
    subroutine Interpolate_cubic1DmonoDerv (nx,f,  fx)
      integer, intent (in)  :: nx
      real,    intent (in)  :: f  (0:nx+1)
      real,    intent (out) :: fx (1:nx  )
    end subroutine Interpolate_cubic1DmonoDerv
  end interface

  interface
    subroutine Interpolate_cubic2Dcoeffs (numberOfSquares, a)
      integer, intent (in)    :: numberOfSquares
      real,    intent (inout) :: a (1:16,*)
    end subroutine Interpolate_cubic2Dcoeffs
  end interface

  interface
    function Interpolate_cubic2DFd1d2 (a,x,y)
      real, intent (in) :: a (1:16)
      real, intent (in) :: x,y
      real :: Interpolate_cubic2DFd1d2 (1:6)
    end function Interpolate_cubic2DFd1d2
  end interface

  interface
    function Interpolate_cubic2DFd1 (a,x,y)
      real, intent (in) :: a (1:16)
      real, intent (in) :: x,y
      real :: Interpolate_cubic2DFd1 (1:3)
    end function Interpolate_cubic2DFd1
  end interface

  interface
    real function Interpolate_cubic2DF (a,x,y)
      real, intent (in) :: a (1:16)
      real, intent (in) :: x,y
    end function Interpolate_cubic2DF
  end interface

  interface
    subroutine Interpolate_cubic2DmonoDerv (nx, ny,        &
                                            f,             &
                                                       fx, &
                                                       fy, &
                                                       fxy )
      integer, intent (in)  :: nx, ny
      real,    intent (in)  :: f   (-1:nx+2 , -1:ny+2)
      real,    intent (out) :: fx  ( 0:nx+1 ,  0:ny+1)
      real,    intent (out) :: fy  ( 0:nx+1 ,  0:ny+1)
      real,    intent (out) :: fxy ( 1:nx   ,  1:ny  )
    end subroutine Interpolate_cubic2DmonoDerv
  end interface

  interface
    subroutine Interpolate_cubic3Dcoeffs (numberOfCubes, a)
      integer, intent (in)    :: numberOfCubes
      real,    intent (inout) :: a (1:64,*)
    end subroutine Interpolate_cubic3Dcoeffs
  end interface

  interface
    function Interpolate_cubic3DFd1d2 (a,x,y,z)
      real, intent (in) :: a (1:64)
      real, intent (in) :: x,y,z
      real :: Interpolate_cubic3DFd1d2 (1:10)
    end function Interpolate_cubic3DFd1d2
  end interface

  interface
    function Interpolate_cubic3DFd1 (a,x,y,z)
      real, intent (in) :: a (1:64)
      real, intent (in) :: x,y,z
      real :: Interpolate_cubic3DFd1 (1:4)
    end function Interpolate_cubic3DFd1
  end interface

  interface
    real function Interpolate_cubic3DF (a,x,y,z)
      real, intent (in) :: a (1:64)
      real, intent (in) :: x,y,z
    end function Interpolate_cubic3DF
  end interface

  interface
    subroutine Interpolate_cubic3DmonoDerv (nx, ny, nz,           &
                                            f,                    &
                                                     fx,fy,fz,    &
                                                     fxy,fxz,fyz, &
                                                     fxyz         )
      integer, intent (in)  :: nx, ny, nz
      real,    intent (in)  :: f    (-2:nx+3 , -2:ny+3 , -2:nz+3)
      real,    intent (out) :: fx   (-1:nx+2 , -1:ny+2 , -1:nz+2)
      real,    intent (out) :: fy   (-1:nx+2 , -1:ny+2 , -1:nz+2)
      real,    intent (out) :: fz   (-1:nx+2 , -1:ny+2 , -1:nz+2)
      real,    intent (out) :: fxy  ( 0:nx+1 ,  0:ny+1 ,  0:nz+1)
      real,    intent (out) :: fxz  ( 0:nx+1 ,  0:ny+1 ,  0:nz+1)
      real,    intent (out) :: fyz  ( 0:nx+1 ,  0:ny+1 ,  0:nz+1)
      real,    intent (out) :: fxyz ( 1:nx   ,  1:ny   ,  1:nz  )
    end subroutine Interpolate_cubic3DmonoDerv
  end interface

  interface
    subroutine Interpolate_finalize ()
    end subroutine Interpolate_finalize
  end interface

  interface
    subroutine Interpolate_init ()
    end subroutine Interpolate_init
  end interface

end Module Interpolate_interface
