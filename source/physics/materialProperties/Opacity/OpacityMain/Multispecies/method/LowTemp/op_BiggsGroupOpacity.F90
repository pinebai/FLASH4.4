!!****if* source/physics/materialProperties/Opacity/OpacityMain/Multispecies/method/LowTemp/op_BiggsGroupOpacity
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
!!  Elower       : the lower energy bound for the group (in eV)
!!  Eupper       : the upper energy bound for the group (in eV)
!!  Opacity      : the Biggs group opacity (in cm^2/g)
!!
!!***
subroutine op_BiggsGroupOpacity (opacityKind,indexElement,Temperature,Elower,Eupper,Opacity)

  use Opacity_data,        ONLY : op_energyDifferenceTolerance
  use Driver_interface,    ONLY : Driver_abortFlash

  use op_lowTempData,      ONLY : op_maxElements,          &
                                  op_maxJmax,              &
                                  op_A1group,              &
                                  op_A2group,              &
                                  op_A3group,              &
                                  op_A4group,              &
                                  op_intLimits,            &
                                  op_elementAij4,          &
                                  op_elementJmax,          &
                                  op_elementPEenergyRange

  use op_numericsData,     ONLY : zero,         &
                                  op_Boltzmann, &
                                  op_keV2erg,   &
                                  op_erg2keV,   &
                                  op_eV2keV

  use op_interface,        ONLY : op_BiggsPlanckGroupIntegrate,  &
                                  op_BiggsRosslndGroupIntegrate

  implicit none

# include "Opacity.h"

  character (len=9), intent (in)  :: opacityKind
  integer,           intent (in)  :: indexElement
  real,              intent (in)  :: Temperature
  real,              intent (in)  :: Elower
  real,              intent (in)  :: Eupper
  real,              intent (out) :: Opacity

  logical :: zeroSectionLower
  logical :: zeroSectionUpper

  integer :: j,n
  integer :: Jlower,Jupper
  integer :: Jmax
  integer :: nSec

  real    :: A1factor, A2factor, A3factor, A4factor
  real    :: BiggsIntegral, Integral
  real    :: Efactor
  real    :: EkeV, ElowerkeV, EupperkeV, ElowestPE, EhighestPE
  real    :: kT
  real    :: rescaleBase10Exponent
!
!
!   ...Check the energy boundaries.
!
!
  if (Eupper <= Elower) then
      call Driver_abortFlash ('[op_BiggsGroupOpacity] ERROR: Bad energy group boundaries')
  end if
!
!
!   ...Convert the energy boundaries to the Biggs units (keV). Abort the calculation, if
!      the lower and upper integration limits are considered to be equal.
!
!
  ElowerkeV = Elower * op_eV2keV
  EupperkeV = Eupper * op_eV2keV

  if (EupperkeV - ElowerkeV < op_energyDifferenceTolerance) then
      call Driver_abortFlash ('[op_BiggsGroupOpacity] ERROR: Integration limits equal')
  end if
!
!
!   ...Decide early, if the energy boundaries are completely outside the PE range either
!      on the low or the high energy side:
!
!                 upper integration energy -> |    | <- lowest Biggs energy
!                      highest Biggs energy-> |    | <- lower integration energy
!
!      If the case, return with opacity equal to zero.
!
!
  Jmax   = op_elementJmax (indexElement)

  ElowestPE  = op_elementPEenergyRange ( LOW,   1,indexElement)
  EhighestPE = op_elementPEenergyRange (HIGH,Jmax,indexElement)

  if (     (EupperkeV - op_energyDifferenceTolerance <= ElowestPE ) &
      .or. (ElowerkeV + op_energyDifferenceTolerance >= EhighestPE) ) then
       Opacity = zero
       return
  end if
!
!
!   ...Identify the Biggs sections of integration.
!
!      First the lower limit. Identify a possible zero section with no contribution
!      to the Biggs integral.
!
!      The following scenarios for the lower limit can happen (in order of testing):
!
!        1) the lower integration limit is below the lowest Biggs table boundary:
!
!                                                  | <- lowest Biggs energy
!                 lower integration energy -> |****
!
!           There is a zone outside the PE range where no radiation is absorbed.
!           This is the zero zone, which is emulated as having zero Biggs expansion
!           coefficients. Note, that ignoring this zone would lead to wrong opacity
!           averaging over the entire energy range. This can be pictured by an extreme
!           example of a large zero zone and a very tiny PE zone, in which case ignorance
!           of the zero zone would lead to large opacity values, but the true situation
!           is a low opacity due to the large zero zone.
!
!        2) the lower integration limit is considered to sit right at the lowest Biggs
!           table boundary (within tolerance):
!
!                                                  | <- lowest Biggs energy
!                      lower integration energy -> |
!
!           There is no zero zone in this case.
!
!        3) the lower integration limit is larger than the lowest Biggs table boundary:
!
!                                                  | <- lowest Biggs energy
!                             lower integration energy -> |
!
!
!
  Jlower = Jmax

  if (ElowerkeV + op_energyDifferenceTolerance < ElowestPE) then
      zeroSectionLower = .true.
      Jlower = 1
  else if (ElowerkeV <= ElowestPE) then
      zeroSectionLower = .false.
      Jlower = 1
  else
      zeroSectionLower = .false.
      do j = 2,Jmax
         EkeV = op_elementPEenergyRange (LOW,j,indexElement)
         if (EkeV > ElowerkeV) then
             Jlower = j - 1
             exit
         else if (EkeV == ElowerkeV) then
             Jlower = j
             exit
         end if
      end do
  end if
!
!
!   ...Next the upper limit.
!
!      The following scenarios for the upper limit can happen (in order of testing):
!
!        1) the upper integration limit is above the highest Biggs table boundary:
!
!                                                  | <- upper integration energy
!                     highest Biggs energy -> |****
!
!           There is a zone outside the PE range where no radiation is absorbed.
!           This is the zero zone, which is emulated as having zero Biggs expansion
!           coefficients. Note, that ignoring this zone would lead to wrong opacity
!           averaging over the entire energy range. This can be pictured by an extreme
!           example of a large zero zone and a very tiny PE zone, in which case ignorance
!           of the zero zone would lead to large opacity values, but the true situation
!           is a low opacity due to the large zero zone.
!
!        2) the upper integration limit is considered to sit right at the highest Biggs
!           table boundary (within tolerance):
!
!                                                  | <- upper integration energy
!                      highest Biggs energy -> |
!
!           There is no zero zone in this case.
!
!        3) the upper integration limit is lower than the highest Biggs table boundary:
!
!                                                  | <- upper integration energy
!                                 highest Biggs energy -> |
!
!
!
  Jupper = 1

  if (EupperkeV - op_energyDifferenceTolerance > EhighestPE) then
      zeroSectionUpper = .true.
      Jupper = Jmax
  else if (EupperkeV >= EhighestPE) then
      zeroSectionUpper = .false.
      Jupper = Jmax
  else
      zeroSectionUpper = .false.
      do j = Jmax-1,1,-1
         EkeV = op_elementPEenergyRange (HIGH,j,indexElement)
         if (EkeV < EupperkeV) then
             Jupper = j + 1
             exit
         else if (EkeV == EupperkeV) then
             Jupper = j
             exit
         end if
      end do
  end if

  if (Jupper < Jlower) then
      call Driver_abortFlash ('[op_BiggsGroupOpacity] ERROR: Biggs limits out of order')
  end if
!
!
!   ...Print section (activate for testing purposes).
!
!
!  write (*,*) ' Jmax             = ',Jmax
!  write (*,*) ' Jlower           = ',Jlower
!  write (*,*) ' Jupper           = ',Jupper
!  write (*,*) ' zeroSectionLower = ',zeroSectionLower
!  write (*,*) ' zeroSectionUpper = ',zeroSectionUpper
!
!
!   ...Group together the needed A1,A2,A3,A4 constants (in opacity units of cm^2/g)
!      and the dimensionless integration limits into sections. The original units
!      of the Ax constants are (cm^2/g)(keV)^x. Add the zero section if necessary.
!
!
  kT       = op_Boltzmann * Temperature           ! in erg
  Efactor  = op_keV2erg / kT                      ! in keV^(-1)
  A1factor = Efactor                              ! in keV^(-1)
  A2factor = A1factor * Efactor                   ! in keV^(-2)
  A3factor = A2factor * Efactor                   ! in keV^(-3)
  A4factor = A3factor * Efactor                   ! in keV^(-4)

  if (zeroSectionLower) then
      nSec = 1
      op_intLimits (1) = ElowerkeV * Efactor
      op_A1group   (1) = zero
      op_A2group   (1) = zero
      op_A3group   (1) = zero
      op_A4group   (1) = zero
  else
      nSec = 0
  end if

  do j = Jlower,Jupper
     nSec = nSec + 1
     EkeV = op_elementPEenergyRange (LOW,j,indexElement)
     op_intLimits (nSec) = max (ElowerkeV, EkeV) * Efactor                   ! dimensionless
     op_A1group   (nSec) = op_elementAij4 (1,j,indexElement) * A1factor      ! in cm^2/g
     op_A2group   (nSec) = op_elementAij4 (2,j,indexElement) * A2factor      ! in cm^2/g
     op_A3group   (nSec) = op_elementAij4 (3,j,indexElement) * A3factor      ! in cm^2/g
     op_A4group   (nSec) = op_elementAij4 (4,j,indexElement) * A4factor      ! in cm^2/g
  end do

  EkeV = op_elementPEenergyRange (HIGH,Jupper,indexElement)
  op_intLimits (nSec+1) = min (EupperkeV, EkeV) * Efactor                    ! dimensionless

  if (zeroSectionUpper) then
      nSec = nSec + 1
      op_intLimits (nSec+1) = EupperkeV * Efactor
      op_A1group   (nSec)   = zero
      op_A2group   (nSec)   = zero
      op_A3group   (nSec)   = zero
      op_A4group   (nSec)   = zero
  end if
!
!
!   ...Print section (activate for testing purposes).
!
!
!  do j = 1,nSec
!     write (*,*) ' j = ',j
!     write (*,*) ' A1 = ',op_A1group (j)
!     write (*,*) ' A2 = ',op_A2group (j)
!     write (*,*) ' A3 = ',op_A3group (j)
!     write (*,*) ' A4 = ',op_A4group (j)
!  end do
!
!  do j = 1,nSec+1
!     write (*,'(A13,2X,ES30.16)') ' int Limit = ',op_intLimits (j)
!  end do
!
!
!   ...Call the corresponding integrator.
!
!
  if (opacityKind == "Planck") then

      call op_BiggsPlanckGroupIntegrate  (nSec,                                &
                                          op_A1group,                          &
                                          op_A2group,                          &
                                          op_A3group,                          &
                                          op_A4group,                          &
                                          op_intLimits,                        &
                                                        rescaleBase10Exponent, &
                                                        BiggsIntegral,         &
                                                        Integral        )

      Opacity = BiggsIntegral / Integral

!      write (*,*) ' Biggs-Planck Integral = ',BiggsIntegral,'x 10 **',rescaleBase10Exponent
!      write (*,*) '       Planck Integral = ',Integral,'x 10 **',rescaleBase10Exponent
!      write (*,*) '               Opacity = ',Opacity

  else if (opacityKind == "Rosseland") then

      call op_BiggsRosslndGroupIntegrate (nSec,                                &
                                          op_A1group,                          &
                                          op_A2group,                          &
                                          op_A3group,                          &
                                          op_A4group,                          &
                                          op_intLimits,                        &
                                                        rescaleBase10Exponent, &
                                                        BiggsIntegral,         &
                                                        Integral               )

      Opacity = Integral / BiggsIntegral

!      write (*,*) ' Biggs-Rosseland Integral = ',BiggsIntegral,'x 10 **',rescaleBase10Exponent
!      write (*,*) '       Rosseland Integral = ',Integral,'x 10 **',rescaleBase10Exponent
!      write (*,*) '                  Opacity = ',Opacity

  else
      call Driver_abortFlash ('[op_BiggsGroupOpacity] ERROR: Unidentified opacity kind')
  end if
!
!
!   ...Ready! 
!
!
  return
end subroutine op_BiggsGroupOpacity
