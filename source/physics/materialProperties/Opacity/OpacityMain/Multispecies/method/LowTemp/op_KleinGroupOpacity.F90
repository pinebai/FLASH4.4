!!****if* source/physics/materialProperties/Opacity/OpacityMain/Multispecies/method/LowTemp/op_KleinGroupOpacity
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

  use Opacity_data,        ONLY : op_atomWeight,           &
                                  op_element2AtomicNumber

  use op_lowTempData,      ONLY : op_maxElements,           &
                                  op_maxJmax,               &
                                  op_A1group,               &
                                  op_A2group,               &
                                  op_A3group,               &
                                  op_A4group,               &
                                  op_intLimits,             &
                                  op_elementAij4,           &
                                  op_elementJmax,           &
                                  op_elementPEenergyRange

  use op_numericsData,     ONLY : op_Boltzmann,             &
                                  op_Avogadro,              &
                                  op_keV2erg,               &
                                  op_erg2keV,               &
                                  op_KleinNishinaPrefactor, &
                                  op_electronRestMassEnergy

  use Driver_interface,    ONLY : Driver_abortFlash

  use op_interface,        ONLY : op_KleinPlanckGroupIntegrate,  &
                                  op_KleinRosslndGroupIntegrate

  implicit none

# include "Opacity.h"

  character (len=9), intent (in)  :: opacityKind
  integer,           intent (in)  :: indexElement
  real,              intent (in)  :: Temperature
  real,              intent (in)  :: Elower
  real,              intent (in)  :: Eupper
  real,              intent (out) :: Opacity

  integer :: Z

  real    :: A,B,L,R,S
  real    :: Efactor
  real    :: ElowerkeV, EupperkeV
  real    :: KleinIntegral, Integral
  real    :: KNintegralFactor
  real    :: kT
  real    :: rescaleBase10Exponent
!
!
!   ...Check the energy boundaries.
!
!
  if (Eupper <= Elower) then
      call Driver_abortFlash ('[op_KleinGroupOpacity] ERROR: Bad energy group boundaries')
  end if
!
!
!   ...Convert the energy boundaries to keV units (temporarily, might change in future).
!
!
!  ElowerkeV = Elower                    ! if Elower in keV
!  EupperkeV = Eupper                    ! if Eupper in keV
  ElowerkeV = Elower * 0.001            ! if Elower in eV
  EupperkeV = Eupper * 0.001            ! if Eupper in eV
!  ElowerkeV = Elower * op_erg2keV       ! if Elower in erg
!  EupperkeV = Eupper * op_erg2keV       ! if Eupper in erg
!
!
!   ...Calculate the Klein-Nishina integral factor and the integral variable factor B.
!
!
  Z                = op_element2AtomicNumber (indexElement)  ! in # of electrons
  A                = op_atomWeight (Z)                       ! in g/(mole of atoms)
  L                = op_KleinNishinaPrefactor                ! in cm^2/(mole of electrons)
!  KNintegralFactor = L / op_Avogadro                         ! in cm^2/electron
  KNintegralFactor = L * (real (Z) / A)                      ! in cm^2/g
  kT               = op_Boltzmann * Temperature              ! in erg
  B                = kT / op_electronRestMassEnergy          ! no units
!
!
!   ...Calculate the dimensionless integration limits.
!
!
  Efactor          = op_keV2erg / kT                         ! in keV^(-1)
  R                = ElowerkeV * Efactor                     ! no units
  S                = EupperkeV * Efactor                     ! no units
!
!
!   ...Call the corresponding integrator.
!
!
  if (opacityKind == "Planck") then

      call op_KleinPlanckGroupIntegrate  (KNintegralFactor,             &
                                          B,                            &
                                          R,S,                          &
                                                 rescaleBase10Exponent, &
                                                 KleinIntegral,         &
                                                 Integral               )

      Opacity = KleinIntegral / Integral

!      write (*,*) ' Klein-Planck Integral = ',KleinIntegral,'x 10 **',rescaleBase10Exponent
!      write (*,*) '       Planck Integral = ',Integral,'x 10 **',rescaleBase10Exponent
!      write (*,*) '               Opacity = ',Opacity

  else if (opacityKind == "Rosseland") then

      call op_KleinRosslndGroupIntegrate (KNintegralFactor,             &
                                          B,                            &
                                          R,S,                          &
                                                 rescaleBase10Exponent, &
                                                 KleinIntegral,         &
                                                 Integral               )

      Opacity = Integral / KleinIntegral

!      write (*,*) ' Klein-Rosseland Integral = ',KleinIntegral,'x 10 **',rescaleBase10Exponent
!      write (*,*) '       Rosseland Integral = ',Integral,'x 10 **',rescaleBase10Exponent
!      write (*,*) '                  Opacity = ',Opacity

  else
      call Driver_abortFlash ('[op_KleinGroupOpacity] ERROR: Unidentified opacity kind')
  end if
!
!
!   ...Ready! 
!
!
  return
end subroutine op_KleinGroupOpacity
