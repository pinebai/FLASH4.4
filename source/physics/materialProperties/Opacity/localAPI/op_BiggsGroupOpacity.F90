!!****if* source/physics/materialProperties/Opacity/localAPI/op_BiggsGroupOpacity
!!
!! NAME
!!
!!  op_BiggsGroupOpacity
!!
!! SYNOPSIS
!!
!!  call op_BiggsGroupOpacity (character (in)  :: opacityKind (len=9),
!!                             integer   (in)  :: indexElement,
!!                             real      (in)  :: Temperature,
!!                             real      (in)  :: Elower,
!!                             real      (in)  :: Eupper,
!!                             real      (out) :: Opacity)
!!
!! DESCRIPTION
!!
!!  This routine constructs the Biggs group (mean) opacity for a particular element.
!!
!! ARGUMENTS
!!
!!  opacityKind  : the kind of Biggs opacity wanted ("Planck" or "Rosseland")
!!  indexElement : the element index handle
!!  Temperature  : the temperature for the Planck radiation law (in K)
!!  Elower       : the lower energy bound for the group (temporarily in keV)
!!  Eupper       : the upper energy bound for the group (temporarily in keV)
!!  Opacity      : the Biggs group opacity (in cm^2/g)
!!
!!***
subroutine op_BiggsGroupOpacity (opacityKind,indexElement,Temperature,Elower,Eupper,Opacity)

  implicit none

  character (len=9), intent (in)  :: opacityKind
  integer,           intent (in)  :: indexElement
  real,              intent (in)  :: Temperature
  real,              intent (in)  :: Elower
  real,              intent (in)  :: Eupper
  real,              intent (out) :: Opacity

  Opacity = 0.0

  return
end subroutine op_BiggsGroupOpacity
