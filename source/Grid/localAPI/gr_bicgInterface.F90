!!****ih* source/Grid/localAPI/gr_bicgInterface
!!
!! NAME
!!  gr_bicgInterface
!!
!! SYNOPSIS
!!  use gr_bicgInterface
!!
!! DESCRIPTION
!! This is the interface module for the BiPCGStab solver implementation.
!!
!! Interfaces for subroutines are defined here.
!! Most of these are meant to be called only within the current (Paramesh) BiPCGStab
!! solver implementation.
!! 
!!***

Module gr_bicgInterface
  implicit none

  interface
  subroutine gr_bicgInit()
    implicit none
  end subroutine gr_bicgInit
  end interface

  interface
  subroutine gr_bicgFinalize()
    implicit none
  end subroutine gr_bicgFinalize
  end interface

  interface
  subroutine gr_bicgInitSlv(bndTypes)
    implicit none
    integer, intent(in) :: bndTypes(6)
  end subroutine gr_bicgInitSlv
  end interface

  interface
  subroutine gr_bicgInitSrc (isrc, isoln, poisfact)
    implicit none
    integer, intent(in) :: isrc, isoln
    real, intent(in)    :: poisfact
  end subroutine gr_bicgInitSrc
  end interface

  interface
  subroutine gr_bicgDotprod (ivar1,ivar2, norm)
    implicit none
    integer, intent(in) :: ivar1,ivar2
    real, intent(inout) :: norm
  end subroutine gr_bicgDotprod
  end interface

  interface
  subroutine gr_bicgNorm (ivar, norm)
    implicit none
    integer, intent(in) :: ivar
    real, intent(inout) :: norm
  end subroutine
  end interface

  interface
  subroutine gr_bicgMultiAx (irhs, ilhs, gcellflg)
    implicit none
    integer, intent(in) :: irhs, ilhs
    logical, intent(in) :: gcellflg
  end subroutine gr_bicgMultiAx
  end interface

  interface
  subroutine gr_bipcgstab (isrc,isoln,poisfact, bcro, bcri, bcvi, bcpi, bcsi, bczi, bcyi, &
                        bc_types,bc_values)
    implicit none
    integer, intent(in)       :: isrc, isoln, bcro, bcri, bcvi, bcpi, bcsi, bczi, bcyi 
    integer, intent(in)       :: bc_types(6)
    real, intent(in), dimension(2,6) :: bc_values
    real          :: poisfact
  end subroutine gr_bipcgstab
  end interface

  interface
  subroutine gr_bicgBndry(ivar, nlayers, gcellflg)
    implicit none
    integer, intent(in) :: ivar, nlayers
    logical, intent(in) :: gcellflg
  end subroutine gr_bicgBndry
  end interface

  interface
  subroutine gr_bicgMapBcType(bcTypeToApply,bcTypeFromGrid,varIndex,gridDataStruct, &
                              axis,face,idest)
    implicit none
    integer, intent(OUT) :: bcTypeToApply
    integer, intent(in) :: bcTypeFromGrid,varIndex,gridDataStruct,axis,face
    integer,intent(IN),OPTIONAL:: idest
  end subroutine gr_bicgMapBcType
  end interface

  interface
  subroutine gr_bicgApplyPrecond(isoln, isrc, bcTypes, bcValues)
    implicit none
    integer, intent(in) :: isoln, isrc
    integer, dimension(6), intent(in) :: bcTypes
    real, dimension(2,6), intent(in)  :: bcValues
  end subroutine gr_bicgApplyPrecond
  end interface

!------------------------------------------------------

end Module gr_bicgInterface
