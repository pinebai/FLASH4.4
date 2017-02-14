!!****if* source/physics/sourceTerms/EnergyDeposition/EnergyDepositionMain/Laser/LaserBeams/ed_beamCrossSectionFunctions
!!
!! NAME
!!
!!  ed_beamCrossSectionFunctions
!!
!! SYNOPSIS
!!
!!  use ed_beamCrossSectionFunctions
!!
!! DESCRIPTION
!!
!!  This module contains all possible beam cross section functions. The
!!  functions return a numerical value depending on coordinate positions
!!  relative to an origin in space. Details are explained below for each
!!  function implemented.
!!
!! ARGUMENTS
!!
!!***


Module ed_beamCrossSectionFunctions

implicit none

contains
!
!
!     ...The main function, which calls one of the others below. This is
!        the wrapper function to the other routines which use it.
!
!
  real function ed_beamCrossSectionWeight (functionType, &
                                           x,y,          &
                                           Cx,Cy,        &
                                           Rx,Ry,        &
                                           Exponent      )
!
!
!     ...The input to this wrapper function consists in specifying the
!        position of the point on a 2D (x,y) area or 1D (x) line.
!
!        The extra variables:
!
!                                Cx,Cy
!                                Rx,Ry
!                                Exponent
!
!        are reals, which are passed for convenient use to be treated as
!        simple factors or as exponents, depending on the function that is
!        being evaluated. If the user feels the need to add more of these
!        reals for more complicated functions, he can do so in form of
!        optional parameters.
!
!        The string variable 'functionType' denotes the type of function
!        used. All other input arguments are optional and the code will check
!        for the particular function selected, if the necessary arguments
!        are present.
!
!        Currently, the set of functions available consists of:
!
!           functionType = 'uniform'            -->  equal weight on entire area
!           functionType = 'gaussian1D'         -->  gaussian weight in 1D
!           functionType = 'gaussian2D'         -->  gaussian weight in 2D
!           functionType = 'gaussianInverse1D'  -->  inverse gaussian weight in 1D
!           functionType = 'gaussianInverse2D'  -->  inverse gaussian weight in 2D
!
!
    use Driver_interface,  ONLY : Driver_abortFlash

    implicit none

    character (len = *), intent (in)            :: functionType
    real,                intent (in), optional  :: x, y
    real,                intent (in), optional  :: Cx, Cy
    real,                intent (in), optional  :: Rx, Ry
    real,                intent (in), optional  :: Exponent

    logical :: argumentsMissing
!
!
!     ...Select the function accordin to its ID.
!
!
    select case (functionType)

    case ('uniform')

      ed_beamCrossSectionWeight = 1.0

    case ('gaussian1D')

      argumentsMissing =     (.not.present (x))        &
                        .or. (.not.present (Cx))       &
                        .or. (.not.present (Rx))       &
                        .or. (.not.present (Exponent))

      if (argumentsMissing) then
          call Driver_abortFlash ("ed_beamCrossSectionWeight: Missing Argument(s)!")
      end if

      ed_beamCrossSectionWeight = ed_gauss1D (x,Cx,Rx,Exponent)

    case ('gaussian2D')

      argumentsMissing =     (.not.present (x))        &
                        .or. (.not.present (y))        &
                        .or. (.not.present (Cx))       &
                        .or. (.not.present (Cy))       &
                        .or. (.not.present (Rx))       &
                        .or. (.not.present (Ry))       &
                        .or. (.not.present (Exponent))

      if (argumentsMissing) then
          call Driver_abortFlash ("ed_beamCrossSectionWeight: Missing Argument(s)!")
      end if

      ed_beamCrossSectionWeight = ed_gauss2D (x,y,Cx,Cy,Rx,Ry,Exponent)

    case ('gaussianInverse1D')

      argumentsMissing =     (.not.present (x))        &
                        .or. (.not.present (Cx))       &
                        .or. (.not.present (Rx))       &
                        .or. (.not.present (Exponent))

      if (argumentsMissing) then
          call Driver_abortFlash ("ed_beamCrossSectionWeight: Missing Argument(s)!")
      end if

      ed_beamCrossSectionWeight = ed_gaussInverse1D (x,Cx,Rx,Exponent)

    case ('gaussianInverse2D')

      argumentsMissing =     (.not.present (x))        &
                        .or. (.not.present (y))        &
                        .or. (.not.present (Cx))       &
                        .or. (.not.present (Cy))       &
                        .or. (.not.present (Rx))       &
                        .or. (.not.present (Ry))       &
                        .or. (.not.present (Exponent))

      if (argumentsMissing) then
          call Driver_abortFlash ("ed_beamCrossSectionWeight: Missing Argument(s)!")
      end if

      ed_beamCrossSectionWeight = ed_gaussInverse2D (x,y,Cx,Cy,Rx,Ry,Exponent)

    case default

      call Driver_abortFlash ("ed_beamCrossSectionWeight: No proper function specified!")

    end select

  end function ed_beamCrossSectionWeight
!
!
!     ...Function(s) for 1D lines.
!
!
!-----------------------------------------------------------------------
  real function ed_gauss1D (x,Cx,Rx,alpha)
!
!
!     ...This function computes the following numerical value:
!
!                             -[({x-Cx}/Rx)^2]^alpha
!                          exp
!
!
    implicit none
    real, intent (in)  :: x
    real, intent (in)  :: Cx
    real, intent (in)  :: Rx
    real, intent (in)  :: alpha
    real :: xf
    xf = (x - Cx) / Rx
    ed_gauss1D = exp ( - (xf * xf) ** alpha )
  end function ed_gauss1D
!-----------------------------------------------------------------------
  real function ed_gaussInverse1D (x,Cx,Rx,alpha)
!
!
!     ...This function computes the following numerical value:
!
!                             +[({x-Cx}/Rx)^2]^alpha
!                          exp
!
!
    implicit none
    real, intent (in)  :: x
    real, intent (in)  :: Cx
    real, intent (in)  :: Rx
    real, intent (in)  :: alpha
    real :: xf
    xf = (x - Cx) / Rx
    ed_gaussInverse1D = exp ( (xf * xf) ** alpha )
  end function ed_gaussInverse1D
!-----------------------------------------------------------------------
!
!
!     ...Function(s) for 2D surfaces.
!
!
!-----------------------------------------------------------------------
  real function ed_gauss2D (x,y,Cx,Cy,Rx,Ry,alpha)
!
!
!     ...This function computes the following numerical value:
!
!                             -[({x-Cx}/Rx)^2 + ({y-Cy}/Ry)^2]^alpha
!                          exp
!
!
    implicit none
    real, intent (in)  :: x, y
    real, intent (in)  :: Cx, Cy
    real, intent (in)  :: Rx, Ry
    real, intent (in)  :: alpha
    real :: xf, yf
    xf = (x - Cx) / Rx
    yf = (y - Cy) / Ry
    ed_gauss2D = exp ( - (xf * xf + yf * yf) ** alpha )
  end function ed_gauss2D
!-----------------------------------------------------------------------
  real function ed_gaussInverse2D (x,y,Cx,Cy,Rx,Ry,alpha)
!
!
!     ...This function computes the following numerical value:
!
!                             +[({x-Cx}/Rx)^2 + ({y-Cy}/Ry)^2]^alpha
!                          exp
!
!
    implicit none
    real, intent (in)  :: x, y
    real, intent (in)  :: Cx, Cy
    real, intent (in)  :: Rx, Ry
    real, intent (in)  :: alpha
    real :: xf, yf
    xf = (x - Cx) / Rx
    yf = (y - Cy) / Ry
    ed_gaussInverse2D = exp ( (xf * xf + yf * yf) ** alpha )
  end function ed_gaussInverse2D
!-----------------------------------------------------------------------

end Module ed_beamCrossSectionFunctions
