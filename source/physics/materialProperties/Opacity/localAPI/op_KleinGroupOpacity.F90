!!****if* source/physics/materialProperties/Opacity/localAPI/op_KleinGroupOpacity
!!
!! NAME
!!
!!  op_KleinGroupOpacity
!!
!! SYNOPSIS
!!
!!  call op_KleinGroupOpacity (character (in)  :: opacityKind (len=9),
!!                             integer   (in)  :: indexElement,
!!                             real      (in)  :: Temperature,
!!                             real      (in)  :: Elower,
!!                             real      (in)  :: Eupper,
!!                             real      (out) :: Opacity)
!!
!! DESCRIPTION
!!
!!  This routine constructs the Klein-Nishina group (mean) opacity for a particular element.
!!
!! ARGUMENTS
!!
!!  opacityKind  : the kind of Klein-Nishina opacity wanted ("Planck" or "Rosseland")
!!  indexElement : the element index handle
!!  Temperature  : the temperature for the Planck radiation law (in K)
!!  Elower       : the lower energy bound for the group (temporarily in eV)
!!  Eupper       : the upper energy bound for the group (temporarily in eV)
!!  Opacity      : the Klein-Nishina group opacity (in cm^2/g)
!!
!!***
subroutine op_KleinGroupOpacity (opacityKind,indexElement,Temperature,Elower,Eupper,Opacity)

  implicit none

  character (len=9), intent (in)  :: opacityKind
  integer,           intent (in)  :: indexElement
  real,              intent (in)  :: Temperature
  real,              intent (in)  :: Elower
  real,              intent (in)  :: Eupper
  real,              intent (out) :: Opacity

  Opacity = 0.0

  return
end subroutine op_KleinGroupOpacity
