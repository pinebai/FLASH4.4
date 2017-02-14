!!****if* source/physics/Eos/EosMain/multiTemp/Gamma/eos_idealGamma
!!
!!
!! NAME
!!
!!  eos_idealGamma
!!
!! SYNOPSIS
!!
!!  call eos_idealGamma(integer(IN) :: mode,
!!                      integer(IN) :: vecLen,
!!                      real(INOUT) :: eosData(vecLen*EOS_NUM),
!!            optional, integer(IN) :: vecBegin,
!!            optional, integer(IN) :: vecEnd,
!!            optional, real(IN)    :: massFrac(vecLen*NSPECIES),
!!      optional,target,logical(IN) :: mask(EOS_VARS+1:EOS_NUM)    )
!!
!! DESCRIPTION
!!
!!
!!  This routine applies the gamma law equation of state to thermodynamic 
!!  quantities at one or more grid points.  The number of points is 
!!  determined by the argument veclen.  Data should be packaged for this 
!!  routine in the 1d array, eosData.  The data in eosData is organized as: 
!!  1:vecLen points contain the first variable, vecLen+1:2*vecLen points 
!!  contain the second variable, and so on. The number and order of
!!  variables in the array is determined by the constants defined in Eos.h.
!!  
!!  The routine takes different quantities as givens depending on the
!!  value of the mode variable: if mode=MODE_DENS_TEMP, density and
!!  temperature are taken as given, and pressure and energy are generated
!!  as output; if mode=MODE_DENS_EI, density and energy are taken as
!!  givens, and pressure and temperature are generated as output.  If
!!  mode=MODE_DENS_PRES, density and pressure are taken as givens, and
!!  energy and temperature are generated as output.
!!  
!!  In addition to pressure, temperature, and internal energy, which are
!!  always thermodynamically consistent after this call, other quantities
!!  such as the various thermodynamic partial derivatives can be
!!  calculated based on the values in the argument, mask.  mask is a
!!  logical array with one entry per quantity, with the order determined
!!  by constants defined in Eos.h (the same as those for the eosData
!!  argument); .true. means return the quantity, .false. means don't.
!!  
!!  This version applies to a single fluid.
!!  
!!  
!! ARGUMENTS 
!! 
!!  mode :    Selects the mode of operation of the Eos unit.
!!             The valid values are MODE_DENS_EI, MODE_DENS_PRES and  
!!             MODE_DENS_TEMP as decribed above.
!!
!!  vecLen   : number of points (cells) for which the eosData array is sized.
!!             If vecBegin and vecEnd are not present, this is also the
!!             number of points (cells) for which EOS computation is to be done.
!!
!!  eosData  : This array is the data structure through which variable values are 
!!             passed in and out of the Eos routine. The arrays is sized as 
!!             EOS_NUM*vecLen. EOS_NUM, and individual input and output
!!             Eos variables are defined in Eos.h. The array is organizes such that
!!             the first 1:vecLen entries represent the first Eos variable, vecLen+1:
!!             2*vecLen represent the second Eos variable and so on. 
!!
!!  vecBegin : Index of first cell in eosData to handle.
!!             Can be used to limit operation to a subrange of cells, untested.
!!             If not present, the default is 1.
!!  vecEnd   : Index of last cell in eosData to handle.
!!             Can be used to limit operation to a subrange of cells, untested.
!!             If not present, the default is vecLen.
!!
!!  massFrac : Contains the mass fractions of the species included in
!!             the simulation. The array is sized as NSPECIES*vecLen.
!!
!!  mask     : Mask is a logical array the size of EOS_DERIVS (number
!!              of partial derivatives that can be computed, defined in
!!              Eos.h), where each index represents a specific partial derivative
!!              that can be calculated by the Eos unit. A .true. value in mask 
!!              results in the corresponding derivative being calculated and 
!!              returned. It should preferably be dimensioned as
!!              mask(EOS_VARS+1:EOS_NUM) in the calling routine 
!!              to exactly match the arguments declaration in Eos Unit.
!!             Note that the indexing of mask does not begin at 1, but rather at one past
!!             the number of variables.
!!
!!             An implementation that does not need derivative quantities should
!!             set the mask equal to .false.
!!
!! PARAMETERS
!!
!!     gamma        :  Ratio of specific heats for the simulated gas
!!     gammaEle     :  Ratio of specific heats for electron component
!!     gammaIon     :  Ratio of specific heats for ion component
!!     eos_singleSpeciesA   :  (Average) atomic mass in amu for simulated gas or plasma
!!     eos_singleSpeciesZ   :  number of electrons per nucleus, i.e., ionization level
!!
!! EXAMPLE
!!
!! --- A single-point at a time example, does not calculate derivatives (based on Cellular Simulation)---
!!
!!  #include "constants.h"   ! for MODE_DENS_TEMP
!!  #include "Flash.h"       ! for NSPECIES
!!  #include "Eos.h"         ! for EOS_VAR order
!!
!!  real  :: temp_zone, rho_zone, ptot, eint, gamma
!!  real, dimension(EOS_NUM)  :: eosData
!!  real, dimension(SPECIES_BEGIN:SPECIES_END) ::  massFraction  
!!  integer, dimension(2,MDIM)                 :: blockRange,blockExtent
!!
!!
!!  massFraction(:) = 1.0e-12        
!!  massFraction(C12_SPEC) = 1.0
!!
!!  .... initiale temp_zone, rho_zone
!!
!!  call Grid_getBlkIndexLimits(blockId,blockRange,blockExtent)
!!  do k = blockRange(LOW,KAXIS), blockRange(HIGH,KAXIS)
!!     do j = blockRange(LOW,JAXIS),blockRange(HIGH,JAXIS)
!!        do i = blockRange(LOW,IAXIS),blockRange(HIGH,IAXIS)
!!
!!           eosData(EOS_TEMP) = temp_zone
!!           eosData(EOS_DENS) = rho_zone
!!
!!           call Eos(MODE_DENS_TEMP,1,eosData,massFraction)
!!
!!           ptot = eosData(EOS_PRES)
!!           eint = eosData(EOS_EINT)
!!           gamma = eosData(EOS_GAMC)
!!
!!           
!!           call Grid_putPointData(blockId,CENTER,TEMP_VAR,EXTERIOR,iPosition,temp_zone)
!!           call Grid_putPointData(blockId,CENTER,DENS_VAR,EXTERIOR,iPosition,rho_zone)
!!           call Grid_putPointData(blockId,CENTER,PRES_VAR,EXTERIOR,iPosition,ptot)
!!           call Grid_putPointData(blockId,CENTER,EINT_VAR,EXTERIOR,iPosition,eint)
!!               if you want ENER_VAR, calculate it from EINT_VAR and kinetic energy
!!           call Grid_putPointData(blockId,CENTER,GAMC_VAR,EXTERIOR,iPosition,gamma)
!!           call Grid_putPointData(blockId,CENTER,GAME_VAR,EXTERIOR,iPosition,(ptot/(etot*sim_rhoAmbient) + 1.0))
!!
!!         enddo  ! end of k loop
!!     enddo     ! end of j loop
!!  enddo        ! end of i loop
!!
!! ------------------ Row at a time example, with derivates (based on Eos_unitTest) --------
!!
!!  #include "constants.h"   ! for MODE_DENS_TEMP
!!  #include "Flash.h"       ! for NSPECIES, EOS_NUM
!!  #include "Eos.h"         ! for EOS_VAR order
!!  integer veclen, isize, jsize, ksize, i,j,k, e
!!  real, dimension(:), allocatable :: eosData
!!  real, dimension(:), allocatable :: massFrac
!!  logical, dimension (EOS_VARS+1:EOS_NUM) :: mask
!!  real, allocatable, dimension(:,:,:,:) :: derivedVariables
!!  integer,dimension(2,MDIM) :: blkLimits,blkLimitsGC
!!
!!   ! in the Eos_unitTest, this loops over all blocks.... here is a snippet from inside
!!     call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)
!!
!!    !  Allocate the necessary arrays for an entire block of data
!!    isize = (blkLimits(HIGH,IAXIS) - blkLimits(LOW,IAXIS) + 1)
!!    jsize = (blkLimits(HIGH,JAXIS) - blkLimits(LOW,JAXIS) + 1)
!!    ksize = (blkLimits(HIGH,KAXIS) - blkLimits(LOW,KAXIS) + 1)
!!    vecLen=isize
!!    allocate(derivedVariables(isize,jsize,ksize,EOS_NUM))
!!    allocate(eosData(vecLen*EOS_NUM))
!!    allocate(massFrac(vecLen*NSPECIES))
!!    mask = .true.
!!
!!    ! indices into the first location for these variables
!!    pres = (EOS_PRES-1)*vecLen
!!    dens = (EOS_DENS-1)*vecLen
!!    temp = (EOS_TEMP-1)*vecLen
!!
!!
!!    call Grid_getBlkPtr(blockID,solnData)
!!    do k = blkLimits(LOW,KAXIS),blkLimits(HIGH,KAXIS)
!!        do j = blkLimits(LOW,JAXIS), blkLimits(HIGH, JAXIS)
!!           do i = 1,vecLen
!!              massFrac((i-1)*NSPECIES+1:i*NSPECIES) = &
!!                   solnData(SPECIES_BEGIN:SPECIES_END,ib+i-1,j,k)
!!           end do
!!
!!           eosData(pres+1:pres+vecLen) =  solnData(PRES_VAR,ib:ie,j,k)
!!           eosData(dens+1:dens+vecLen) =  solnData(DENS_VAR,ib:ie,j,k)
!!           ! Eos Helmholtz needs a good initial estimate of temperature no matter what the mode
!!           eosData(temp+1:temp+vecLen) =  solnData(TEMP_VAR,ib:ie,j,k)
!!
!!           call Eos(MODE_DENS_PRES,vecLen,eosData,massFrac,mask)
!!
!!           do e=EOS_VARS+1,EOS_NUM
!!              m = (e-1)*vecLen
!!              derivedVariables(1:vecLen,j-NGUARD,k-NGUARD,e) =  eosData(m+1:m+vecLen)
!!           end do
!!        end do
!!     end do
!!
!! NOTES
!!
!!  NSPECIES is defined in Flash.h.
!!
!!  EOS_VARS and EOS_NUM  are defined in Eos.h.
!!  Calling funtions should included Eos.h, in order to get the definitions of
!!  Eos-specific constants to be able to populate the eosData and mask arrays.
!!  
!!  MODE_DENS_TEMP, MODE_DENS_EI, MODE_DENS_PRES, etc. are defined in constants.h.
!!
!!  All routines calling this routine should include a 
!!  use Eos_interface 
!!  statement, preferable with "ONLY" attribute e.g.
!!  use Eos_interface, ONLY:  Eos
!!
!!  For Gamma and Multigamma routines, the entropy and entropy derivatives 
!!  calculations have not been confirmed to be correct.  Use with caution.
!!  In THIS implementation, we do not even TRY to return combined entropies.
!!  Hopefully that will change at some point.
!!
!! SEE ALSO
!! 
!!  Eos.h    defines the variables used.
!!  Eos_wrapped  sets up the data structure.
!!
!!***
!#ifdef DEBUG_ALL
#define DEBUG_EOS
!#endif

subroutine eos_idealGamma(mode, vecLen, eosData, vecBegin,vecEnd, massFrac, mask, diagFlag)

!==============================================================================
  use Eos_data, ONLY : eos_gasConstant, &
       eos_singleSpeciesA, eos_singleSpeciesZ, &
       eos_smallT, eos_smallEion, eos_smallEele, eos_smallRho, eos_tol, eos_maxNewton, &
       eos_forceConstantInput, &
       eos_meshMe, &
       eos_combinedTempRule, eos_entrEleScaleChoice
  use eos_helmConstData, ONLY: eos_ao3, eos_kerg, eos_avo, eos_kergavo, eos_h,eos_hbar
  use eos_idealGammaData, ONLY: eos_gamma, eos_gammam1, &
       eos_gammaIon,eos_gammaEle, &
       eos_gammam1Ion, eos_gammam1Ele, eos_gammam1Rad, &
       eos_eMass, eos_eMassInUAmu
  use Driver_interface, ONLY : Driver_abortFlash
  use Hydro_interface, ONLY : Hydro_recalibrateEints

  use eos_mtInterface, ONLY : eos_byTemp

  use eos_vecData, ONLY:  tempRow,&
       denRow, abarRow, zbarRow, &
       gamcRow, ptotRow, stotRow, eCompRow, pCompRow, &
       detRow, dptRow, dpdRow, dedRow, pelRow, neRow, etaRow, cvRow, cpRow, &
       tempIonRow,tempEleRow

#include "Eos_components.h"
  implicit none

#include "constants.h"
#include "Eos.h"
#include "Flash.h"

  !     Arguments
  integer, INTENT(in) :: mode, vecLen
  real,INTENT(inout), dimension(EOS_NUM*vecLen) :: eosData 
  integer,optional,INTENT(in) :: vecBegin,vecEnd
  real, optional, INTENT(in),dimension(NSPECIES*vecLen)    :: massFrac
  logical, optional, INTENT(in),target,dimension(EOS_VARS+1:EOS_NUM) :: mask
  integer, optional, INTENT(out)    :: diagFlag

  ! This is the variable that is used internally -- set to false unless mask comes in
  logical,save,target, dimension(EOS_VARS+1:EOS_NUM) :: maskInternal
  data maskInternal / EOS_DERIVS * .FALSE. /
  logical,pointer, dimension(:) :: maskPtr

#ifdef EOS_FORCE_2T
  integer,parameter :: overallCMask(EOSCOMP_NUM_COMPONENTS) = (/1,1,0/)
#else
  integer,parameter :: overallCMask(EOSCOMP_NUM_COMPONENTS) = (/1,1,1/)
#endif
  real,parameter :: eos_gammaRad = (4./3.)

  integer :: componentMask(EOSCOMP_NUM_COMPONENTS)

  real ::  ggprod, ggprodinv, gam1inv
  real :: ggprodIon, ggprodEle, ggprodRad
  real :: ggprodinvIon, ggprodinvEle, ggprodinvRad
  real :: gam1invIon, gam1invEle, gam1invRad
  real :: ionFrac, eleFrac, radFrac
  real :: dynamicZ, relA, Zinv, Zp, ZpInv
  real :: Ye
  real :: kBoltz
  real :: tempRef,tempRefInv
  real :: densRef,densRefInv
  integer :: dens, temp, pres, eint, abar, zbar
  integer :: entr, dst, dsd
  integer :: dpt, dpd, det, ded, c_v, c_p, gamc, pel, ne, eta
  integer :: tempIon,tempEle,tempRad                                                
  integer :: eintIon,eintEle,eintRad                                                
  integer :: presIon,presEle,presRad
  integer :: entrEle,entrRad
  integer :: i, ilo,ihi, rowLen
  integer :: seleVariant
  logical :: arrayBoundHackMode
  logical :: useNRForAddtlCombined !use vectors set by NR call for addtl. combined components
  !addtl. comp. = cp, cv, ...
  logical :: useNRForGamcCombined  !use vectors set by NR call for gamc combined components
  logical :: doUncoupledForAddtl ! assume components thermally uncoupled for addtl. out variables
  logical :: doUncoupledForGamc  ! assume components thermally uncoupled for gamc out variables

  maskPtr => maskInternal
  if (present(mask)) then
     maskPtr => mask
  end if
 

  if (lbound(maskPtr,1) == 1) then
     !gfortran does not set lower bound of pointer to EOS_VARS+1.
     !(I think gfortran is doing the correct thing).
     arrayBoundHackMode = .true.

     nullify(maskPtr)
     allocate(maskPtr(EOS_VARS+1:EOS_NUM))
     if (present(mask)) then
        maskPtr(EOS_VARS+1:EOS_NUM) = mask(EOS_VARS+1:EOS_NUM)
     else
        maskPtr(EOS_VARS+1:EOS_NUM) = maskInternal(EOS_VARS+1:EOS_NUM)
     end if
  else
     arrayBoundHackMode = .false.
  end if

  if (present(diagFlag)) diagFlag = 0

  useNRForAddtlCombined = .FALSE.; useNRForGamcCombined = .FALSE.
  doUncoupledForAddtl = .FALSE.; doUncoupledForGamc = .FALSE.

  kBoltz = eos_kerg
  seleVariant = eos_entrEleScaleChoice        !Use 4,5,6,7,8; other values experimental/historical - KW

  densRef = eos_smallRho; densRefInv = 1.0/densRef
  tempRef = 1.0; tempRef = eos_smallT; tempRefInv = 1.0/tempRef

  ggprod = eos_gammam1 * eos_gasConstant


  !We split up the total fluid density in the following way:

  dynamicZ = eos_singleSpeciesZ
  Zp = dynamicZ + 1
  Zinv = 1.0/dynamicZ; ZpInv = 1.0/Zp
!  nucZ = eos_singleSpeciesZ     !unused
  relA = eos_singleSpeciesA
  eleFrac = dynamicZ*eos_eMassInUAmu / relA !fraction of the mass in cell that is in electrons
  ionFrac = 1.0 - eleFrac                   !fraction of the mass in cell that is in ions

  radFrac = 0.0
  if ( ((ionFrac + eleFrac + radFrac) < 0.999999) .or. &
       ((ionFrac + eleFrac + radFrac) > 1.000001) ) then
     call Driver_abortFlash("Fractional division of fluid variables /= 1")
  end if
  

  select case (mode)
  case(MODE_DENS_TEMP,MODE_DENS_TEMP_COMP,MODE_DENS_TEMP_EQUI, &
       MODE_DENS_TEMP_ALL,MODE_DENS_TEMP_GATHER)
     ggprodIon = eos_gammam1Ion * eos_gasConstant
     ggprodEle = eos_gammam1Ele * eos_gasConstant
     ggprodRad = eos_gammam1Rad * radFrac * eos_gasConstant

  ! density, internal energies (or electron entropy) taken as input:
  case(MODE_DENS_EI, MODE_DENS_EI_SCATTER, &
       MODE_DENS_EI_ION, MODE_DENS_EI_ELE, &
       MODE_DENS_EI_COMP, MODE_DENS_EI_GATHER, &
       MODE_DENS_EI_RECAL_GATHER, &
       MODE_DENS_EI_ALL, &
       MODE_DENS_EI_MAT_GATHER,MODE_DENS_EI_MAT_EQUI, &
       MODE_DENS_EI_SELE_GATHER)
     ggprodIon = eos_gammam1Ion * eos_gasConstant
     ggprodEle = eos_gammam1Ele * eos_gasConstant
     ggprodRad = eos_gammam1Rad * eos_gasConstant

     ggprodinv = 1. / ggprod
     gam1inv   = 1. / eos_gammam1
     ggprodinvIon = 1. / ggprodIon
     gam1invIon   = 1. / eos_gammam1Ion
     ggprodinvEle = 1. / ggprodEle
     gam1invEle   = 1. / eos_gammam1Ele
     ggprodinvRad = 1. / ggprodRad
     gam1invRad   = 1. / eos_gammam1Rad


  ! density, pressures taken as input:
  case (MODE_DENS_PRES)
  case default
     call Driver_abortFlash("[eos_idealGamma] Unrecognized input mode given to Eos")
  end select




!============================================================================

  pres = (EOS_PRES-1)*vecLen
  dens = (EOS_DENS-1)*vecLen
  temp = (EOS_TEMP-1)*vecLen
  eint = (EOS_EINT-1)*vecLen
  gamc = (EOS_GAMC-1)*vecLen
  abar = (EOS_ABAR-1)*vecLen
  zbar = (EOS_ZBAR-1)*vecLen
  entr = (EOS_ENTR-1)*vecLen
  tempIon = (EOS_TEMPION-1)*vecLen
  tempEle = (EOS_TEMPELE-1)*vecLen
  tempRad = (EOS_TEMPRAD-1)*vecLen
  eintIon = (EOS_EINTION-1)*vecLen
  eintEle = (EOS_EINTELE-1)*vecLen
  eintRad = (EOS_EINTRAD-1)*vecLen
  presIon = (EOS_PRESION-1)*vecLen
  presEle = (EOS_PRESELE-1)*vecLen
  presRad = (EOS_PRESRAD-1)*vecLen
  entrEle = (EOS_ENTRELE-1)*veclen
  entrRad = (EOS_ENTRRAD-1)*veclen
  

  if (present(vecBegin)) then
     ilo = vecBegin
  else
     ilo = 1
  end if
  if (present(vecEnd)) then
     ihi = vecEnd
  else
     ihi = vecLen
  end if
  rowLen = ihi - ilo + 1
#ifdef DEBUG_EOS
  if (ilo < 1 .OR. ilo > vecLen) then
     print*,'[eos_idealGamma] ilo is',ilo
     call Driver_abortFlash("[eos_idealGamma] invalid ilo")
  end if
  if (ihi < 1 .OR. ihi > vecLen) then
     print*,'[eos_idealGamma] ihi is',ihi
     call Driver_abortFlash("[eos_idealGamma] invalid ihi")
  end if
  if (rowLen < 0 .OR. rowLen > vecLen) then
     print*,'[eos_idealGamma] rowLen is',rowLen
     call Driver_abortFlash("[eos_idealGamma] invalid rowLen")
  end if
#endif
  if (rowLen == 0) then
     print*,'[eos_idealGamma] rowLen is 0.'
  end if

  !! For allocatable arrays, set them up now.
#ifndef FIXEDBLOCKSIZE
  call eos_vecAlloc(rowLen)
#endif


!HISTORICAL-NOTE  But for NSPECIES = 1, abar is NOT just 1.
!HISTORICAL-NOTE The Flash2 routines initialize a fake "fluid"
!HISTORICAL-NOTE with properties A=1, Z=1, and Eb=1.  Therefore flash2 can
!HISTORICAL-NOTE  always call the equivalent of Multispecies_getSumInV and it
!HISTORICAL-NOTE returns a default number.   In FLASH3, the decision has been made
!HISTORICAL-NOTE that Eos/Gamma cannot be run with more than one fluid.

  do i = ilo,ihi
     eosData(gamc+i) = eos_gamma
     eosData(abar+i) = eos_singleSpeciesA
     eosData(zbar+i) = eos_singleSpeciesZ
     if (mode==MODE_DENS_EI_RECAL_GATHER) then
#ifdef EOS_FORCE_2T
        call Hydro_recalibrateEints(eosData(eint+i), &
             eosData(eintIon+i),eosData(eintEle+i))
#else
        call Hydro_recalibrateEints(eosData(eint+i), &
             eosData(eintIon+i),eosData(eintEle+i),eosData(eintRad+i))
#endif
     end if
  end do
  Ye = eos_singleSpeciesZ / eos_singleSpeciesA
!!$  print 9999,'Ye = ',Ye
!!$  9999 format(A7,1PG25.19)

#define Ye_ARRAYREFERENCE Ye

  ! density, temperature taken as input
  select case (mode)
  case (MODE_DENS_TEMP)
     eosData(pres+ilo:pres+ihi) = eos_gasConstant*eosData(dens+ilo:dens+ihi) * &
          eosData(temp+ilo:temp+ihi) / eosData(abar+ilo:abar+ihi) * Zp
     eosData(eint+ilo:eint+ihi) = ggprod * eosData(temp+ilo:temp+ihi) &
          / eosData(abar+ilo:abar+ihi) * Zp

     !Currently causes an FPE  (entr not used anywhere so just set to 0 for now)
     !eosData(entr+ilo:entr+ihi) = (eosData(pres+ilo:pres+ihi)/eosData(dens+ilo:dens+ihi) +  &
     !     &  eosData(eint+ilo:eint+ihi))/eosData(temp+ilo:temp+ihi)
     eosData(entr+ilo:entr+ihi) = 0.0


  case (MODE_DENS_TEMP_COMP,MODE_DENS_TEMP_GATHER)
     useNRForAddtlCombined = .FALSE.
     doUncoupledForAddtl = .TRUE.
     eosData(presIon+ilo:presIon+ihi) = eos_gasConstant * &
          eosData(dens+ilo:dens+ihi) * &
          eosData(tempIon+ilo:tempIon+ihi) / eosData(abar+ilo:abar+ihi)
     eosData(presEle+ilo:presEle+ihi) = eos_gasConstant * &
          eosData(dens+ilo:dens+ihi) * &
          eosData(tempEle+ilo:tempEle+ihi) / eosData(abar+ilo:abar+ihi) * dynamicZ

     eosData(eintIon+ilo:eintIon+ihi) = ggprodIon * eosData(tempIon+ilo:tempIon+ihi) &
          / eosData(abar+ilo:abar+ihi)
     eosData(eintEle+ilo:eintEle+ihi) = ggprodEle * eosData(tempEle+ilo:tempEle+ihi) &
          / eosData(abar+ilo:abar+ihi) * dynamicZ

     eosData(presRad+ilo:presRad+ihi)    = eos_ao3 * eosData(tempRad+ilo:tempRad+ihi)**4

     eosData(eintRad+ilo:eintRad+ihi)    = 3.0e0 * eosData(presRad+ilo:presRad+ihi) / eosData(dens+ilo:dens+ihi)

     eosData(entr+ilo:entr+ihi) = 0.0
     if (mode==MODE_DENS_TEMP_GATHER) then
!!$        if (ANY(maskPtr((/EOS_PRESELE,EOS_PRESION,EOS_PRESRAD/)) .EQV. .FALSE.)) then
!!$           call Driver_abortFlash("[eos_idealGamma] cannot calculate MODE_DENS_TEMP_GATHER without component pressure masks.&
!!$                & Set mask appropriately.")
!!$        end if
        eosData(eint+ilo:eint+ihi) = ( &
             eosData(eintIon+ilo:eintIon+ihi)+ &
             eosData(eintEle+ilo:eintEle+ihi)+ &
             eosData(eintRad+ilo:eintRad+ihi))

        if(eos_combinedTempRule==0) then
           call eos_newtonRaphson(vecLen, MODE_DENS_EI, .TRUE., maskPtr,cMask=overallCMask)
        endif

        ! Next lines moved from end of eos_newtonRaphson; pres changed - KW 2010-12-13
        eosData(pres+ilo:pres+ihi) = ( &
             eosData(presIon+ilo:presIon+ihi)+ &
             eosData(presEle+ilo:presEle+ihi)+ &
             eosData(presRad+ilo:presRad+ihi))
#ifdef EOS_FORCE_2T
        eosData(gamc+ilo:gamc+ihi)=eos_gammaEle
#else
        call setCombinedGamc(gamc)
!        eosData(gamc+ilo:gamc+ihi)=gamcRow(1:rowLen) !replace with expression based on component states!
#endif
        call setCombinedTemp()
        eosData(entr+ilo:entr+ihi)=stotRow(1:rowLen) !(replace with expression based on component states)

     end if

    !! Note that we make here a choice for the arbitrary additive constant in the electron entropy 
     if(maskPtr(EOS_ENTRELE)) &
          call setEleEntropy(entrele,tempEle)
     if(maskPtr(EOS_ENTRRAD)) &
          call setRadEntropy(entrrad,tempRad,eintRad)
     
        
  case (MODE_DENS_TEMP_EQUI)
     !Equalize temperature.  Each component uses same temperature.
     eosData(tempIon+ilo:tempIon+ihi) = eosData(temp+ilo:temp+ihi)
     eosData(tempEle+ilo:tempEle+ihi) = eosData(temp+ilo:temp+ihi)
     eosData(tempRad+ilo:tempRad+ihi) = eosData(temp+ilo:temp+ihi)

     eosData(pres+ilo:pres+ihi) = eos_gasConstant*eosData(dens+ilo:dens+ihi) * &
          eosData(temp+ilo:temp+ihi) / eosData(abar+ilo:abar+ihi) * Zp
     eosData(presIon+ilo:presIon+ihi) = eos_gasConstant * &
          eosData(dens+ilo:dens+ihi) * &
          eosData(temp+ilo:temp+ihi) / eosData(abar+ilo:abar+ihi)
     eosData(presEle+ilo:presEle+ihi) = eos_gasConstant * &
          eosData(dens+ilo:dens+ihi) * &
          eosData(temp+ilo:temp+ihi) / eosData(abar+ilo:abar+ihi) * dynamicZ
     
     eosData(eint+ilo:eint+ihi) = ggprod * eosData(temp+ilo:temp+ihi) &
          / eosData(abar+ilo:abar+ihi) * Zp
     eosData(eintIon+ilo:eintIon+ihi) = ggprodIon * eosData(temp+ilo:temp+ihi) &
          / eosData(abar+ilo:abar+ihi)
     eosData(eintEle+ilo:eintEle+ihi) = ggprodEle * eosData(temp+ilo:temp+ihi) &
          / eosData(abar+ilo:abar+ihi) * dynamicZ

     eosData(presRad+ilo:presRad+ihi)    = eos_ao3 * eosData(tempRad+ilo:tempRad+ihi)**4

     eosData(eintRad+ilo:eintRad+ihi)    = 3.0e0 * eosData(presRad+ilo:presRad+ihi) / eosData(dens+ilo:dens+ihi)

     eosData(pres+ilo:pres+ihi) = eosData(pres+ilo:pres+ihi) + eosData(presRad+ilo:presRad+ihi)
     eosData(eint+ilo:eint+ihi) = eosData(eint+ilo:eint+ihi) + eosData(eintRad+ilo:eintRad+ihi)

     eosData(entr+ilo:entr+ihi) = 0.0

    !! Note that we make here a choice for the arbitrary additive constant in the electron entropy 
     if(maskPtr(EOS_ENTRELE)) &
          call setEleEntropy(entrele,temp)
     if(maskPtr(EOS_ENTRRAD)) &
          call setRadEntropy(entrrad,temp,eintRad)


  case (MODE_DENS_TEMP_ALL)
     eosData(pres+ilo:pres+ihi) = eos_gasConstant*eosData(dens+ilo:dens+ihi) * &
          eosData(temp+ilo:temp+ihi) / eosData(abar+ilo:abar+ihi) * Zp

     eosData(presIon+ilo:presIon+ihi) = eos_gasConstant * &
          eosData(dens+ilo:dens+ihi) * &
          eosData(tempIon+ilo:tempIon+ihi) / eosData(abar+ilo:abar+ihi)
     eosData(presEle+ilo:presEle+ihi) = eos_gasConstant * &
          eosData(dens+ilo:dens+ihi) * &
          eosData(tempEle+ilo:tempEle+ihi) / eosData(abar+ilo:abar+ihi) * dynamicZ

     eosData(eint+ilo:eint+ihi) = ggprod * eosData(temp+ilo:temp+ihi) &
          / eosData(abar+ilo:abar+ihi) * Zp
     eosData(eintIon+ilo:eintIon+ihi) = ggprodIon * eosData(tempIon+ilo:tempIon+ihi) &
          / eosData(abar+ilo:abar+ihi)
     eosData(eintEle+ilo:eintEle+ihi) = ggprodEle * eosData(tempEle+ilo:tempEle+ihi) &
          / eosData(abar+ilo:abar+ihi) * dynamicZ

     eosData(presRad+ilo:presRad+ihi)    = eos_ao3 * eosData(tempRad+ilo:tempRad+ihi)**4

     eosData(eintRad+ilo:eintRad+ihi)    = 3.0e0 * eosData(presRad+ilo:presRad+ihi) / eosData(dens+ilo:dens+ihi)

     eosData(entr+ilo:entr+ihi) = 0.0

    !! Note that we make here a choice for the arbitrary additive constant in the electron entropy 
     if(maskPtr(EOS_ENTRELE)) &
          call setEleEntropy(entrele,tempEle)
     if(maskPtr(EOS_ENTRRAD)) &
          call setRadEntropy(entrrad,tempRad,eintRad)


  ! density, internal energy taken as input

  case (MODE_DENS_EI, MODE_DENS_EI_SCATTER)
     useNRForAddtlCombined = .TRUE.
     doUncoupledForAddtl = .FALSE.
     eosData(pres+ilo:pres+ihi) = eosData(dens+ilo:dens+ihi) * &
          eosData(eint+ilo:eint+ihi) * gam1inv                        !overridden below!
     eosData(temp+ilo:temp+ihi) = eosData(eint+ilo:eint+ihi) * ggprodinv * &
          eosData(abar+ilo:abar+ihi) * ZpInv                          !overridden below!
     eosData(entr+ilo:entr+ihi) = 0.0                                 !overridden below!

     call eos_newtonRaphson(vecLen, MODE_DENS_EI, .TRUE., maskPtr,cMask=overallCMask)
     
     ! Next lines moved from end of eos_newtonRaphson
     eosData(pres+ilo:pres+ihi)=ptotRow(1:rowLen)
     eosData(gamc+ilo:gamc+ihi)=gamcRow(1:rowLen)
     eosData(entr+ilo:entr+ihi)=stotRow(1:rowLen)

     if (mode == MODE_DENS_EI_SCATTER) then
        eosData(tempIon+ilo:tempIon+ihi)=tempRow(1:rowLen)
        eosData(tempEle+ilo:tempEle+ihi)=tempRow(1:rowLen)
        eosData(tempRad+ilo:tempRad+ihi)=tempRow(1:rowLen)
        if(maskPtr(EOS_PRESION)) eosData(presIon+ilo:presIon+ihi)=pCompRow(EOSCOMP_ION,1:rowLen)
        if(maskPtr(EOS_PRESELE)) eosData(presEle+ilo:presEle+ihi)=pCompRow(EOSCOMP_ELE,1:rowLen)
        if(maskPtr(EOS_PRESRAD)) eosData(presRad+ilo:presRad+ihi)=pCompRow(EOSCOMP_RAD,1:rowLen)
        if(maskPtr(EOS_EINTION)) eosData(eintIon+ilo:eintIon+ihi)=eCompRow(EOSCOMP_ION,1:rowLen)
        if(maskPtr(EOS_EINTELE)) eosData(eintEle+ilo:eintEle+ihi)=eCompRow(EOSCOMP_ELE,1:rowLen)
        if(maskPtr(EOS_EINTRAD)) eosData(eintRad+ilo:eintRad+ihi)=eCompRow(EOSCOMP_RAD,1:rowLen)
     end if


  case(MODE_DENS_EI_ION,MODE_DENS_EI_ELE,MODE_DENS_EI_RAD)

     componentMask(:) = 0
     select case (mode)
     case(MODE_DENS_EI_RAD)
        componentMask(EOSCOMP_RAD) = 1
     case(MODE_DENS_EI_ELE)
        componentMask(EOSCOMP_ELE) = 1
     case(MODE_DENS_EI_ION)
        componentMask(EOSCOMP_ION) = 1
     end select

     select case (mode)
     case(MODE_DENS_EI_RAD)
        eosData(tempRad+ilo:tempRad+ihi) = sqrt(sqrt( eosData(presRad+ilo:presRad+ihi) / eos_ao3 ))
        if(maskPtr(EOS_PRESRAD)) &
             eosData(presRad+ilo:presRad+ihi) = eosData(eintRad+ilo:eintRad+ihi) * eosData(dens+ilo:dens+ihi) / 3.0

     case(MODE_DENS_EI_ELE)
        eosData(tempEle+ilo:tempEle+ihi) = eosData(eintEle+ilo:eintEle+ihi) * ggprodinvEle * &
             eosData(abar+ilo:abar+ihi) * Zinv
        if(any(maskPtr(EOS_VARS+1:EOS_NUM))) then
           useNRForAddtlCombined = .TRUE.; useNRForGamcCombined = .TRUE.
           tempRow(1:rowLen)    = eosData(tempEle+ilo:tempEle+ihi)
           denRow(1:rowLen)     = eosData(dens+ilo:dens+ihi)
           abarRow(1:rowLen)    = eosData(abar+ilo:abar+ihi)
           zbarRow(1:rowLen)    = eosData(zbar+ilo:zbar+ihi)
           tempEleRow(1:rowLen) = tempRow(1:rowLen)
           tempIonRow(1:rowLen) = tempRow(1:rowLen)
           call eos_byTemp(1,rowLen,maskPtr,componentMask=componentMask,eleFrac=eleFrac,ggprodEle=ggprodEle)
           if(maskPtr(EOS_PRESELE)) eosData(presEle+ilo:presEle+ihi) = pCompRow(EOSCOMP_ELE,1:rowLen)
           if(maskPtr(EOS_EINTELE) .AND. .NOT.eos_forceConstantInput) &
                                    eosData(eintEle+ilo:eintEle+ihi) = eCompRow(EOSCOMP_ELE,1:rowLen)
        else if(maskPtr(EOS_PRESELE)) then
           eosData(presEle+ilo:presEle+ihi) = eosData(dens+ilo:dens+ihi) * &
                                              eosData(eintEle+ilo:eintEle+ihi) * gam1invEle
        end if

     case(MODE_DENS_EI_ION)
        eosData(tempIon+ilo:tempIon+ihi) = eosData(eintIon+ilo:eintIon+ihi) * ggprodinvIon * &
             eosData(abar+ilo:abar+ihi)
        if(any(maskPtr(EOS_VARS+1:EOS_NUM))) then
           useNRForAddtlCombined = .TRUE.; useNRForGamcCombined = .TRUE.
           tempRow(1:rowLen)    = eosData(tempIon+ilo:tempIon+ihi)
           denRow(1:rowLen)     = eosData(dens+ilo:dens+ihi)
           abarRow(1:rowLen)    = eosData(abar+ilo:abar+ihi)
           zbarRow(1:rowLen)    = eosData(zbar+ilo:zbar+ihi)
           tempEleRow(1:rowLen) = tempRow(1:rowLen)
           tempIonRow(1:rowLen) = tempRow(1:rowLen)
           call eos_byTemp(1,rowLen,maskPtr,componentMask=componentMask,eleFrac=eleFrac,ggprodEle=ggprodEle)
           if(maskPtr(EOS_PRESION)) eosData(presIon+ilo:presIon+ihi) = pCompRow(EOSCOMP_ION,1:rowLen)
           if(maskPtr(EOS_EINTION) .AND. .NOT.eos_forceConstantInput) &
                                    eosData(eintIon+ilo:eintIon+ihi) = eCompRow(EOSCOMP_ION,1:rowLen)
        else if(maskPtr(EOS_PRESION)) then
           eosData(presIon+ilo:presIon+ihi) = eosData(dens+ilo:dens+ihi) * &
                                              eosData(eintIon+ilo:eintIon+ihi) * gam1invIon
        end if
     end select


  case (MODE_DENS_EI_COMP, MODE_DENS_EI_GATHER, MODE_DENS_EI_RECAL_GATHER, MODE_DENS_EI_ALL)
     if (ANY(abs(eosData(eint+ilo:eint+ihi) - ( &
             eosData(eintIon+ilo:eintIon+ihi)+ &
             eosData(eintEle+ilo:eintEle+ihi)+ &
             eosData(eintRad+ilo:eintRad+ihi))) .GE. 1.e-10*abs(eosData(eint+ilo:eint+ihi)))) then
#ifdef DEBUG_MAR2012
        if (present(diagFlag)) then
           diagflag = 4
        end if
        print*,'EINT1',eosData(eint+ilo:eint+ihi)
        print*,'ESUM1',(eosData(eintIon+ilo:eintIon+ihi)+eosData(eintEle+ilo:eintEle+ihi))
        print*,'EION1',eosData(eintIon+ilo:eintIon+ihi)
        print*,'EELE1',eosData(eintEle+ilo:eintEle+ihi)
#endif
     end if
     useNRForAddtlCombined = .FALSE.
     doUncoupledForAddtl = .TRUE. !or FALSE for _COMP?
     eosData(presIon+ilo:presIon+ihi) =      eosData(dens+ilo:dens+ihi)  * &
          eosData(eintIon+ilo:eintIon+ihi) * gam1invIon
     eosData(presEle+ilo:presEle+ihi) =      eosData(dens+ilo:dens+ihi)  * &
          eosData(eintEle+ilo:eintEle+ihi) * gam1invEle

     eosData(tempIon+ilo:tempIon+ihi) = eosData(eintIon+ilo:eintIon+ihi) * ggprodinvIon * &
          eosData(abar+ilo:abar+ihi)
     eosData(tempEle+ilo:tempEle+ihi) = eosData(eintEle+ilo:eintEle+ihi) * ggprodinvEle * &
          eosData(abar+ilo:abar+ihi) * Zinv

     eosData(presRad+ilo:presRad+ihi) = eosData(eintRad+ilo:eintRad+ihi) * eosData(dens+ilo:dens+ihi) / 3.0
     eosData(tempRad+ilo:tempRad+ihi) = sqrt(sqrt( eosData(presRad+ilo:presRad+ihi) / eos_ao3 ))

     !Currently causes an FPE  (entr not used anywhere so just set to 0 for now)
     !eosData(entr+ilo:entr+ihi) = ( &
     !     (eosData(presIon+ilo:presIon+ihi)/eosData(dens+ilo:dens+ihi) +  &
     !     &  eosData(eintIon+ilo:eintIon+ihi))/eosData(tempIon+ilo:tempIon+ihi) &
     !     ) + ( &
     !     (eosData(presEle+ilo:presEle+ihi)/eosData(dens+ilo:dens+ihi) +  &
     !     &  eosData(eintEle+ilo:eintEle+ihi))/eosData(tempEle+ilo:tempEle+ihi) &
     !     ) + ( &
     !     (eosData(presRad+ilo:presRad+ihi)/eosData(dens+ilo:dens+ihi) +  &
     !     &  eosData(eintRad+ilo:eintRad+ihi))/eosData(tempRad+ilo:tempRad+ihi) &
     !     )
     eosData(entr+ilo:entr+ihi) = 0.0


     if (mode==MODE_DENS_EI_ALL) then
        useNRForAddtlCombined = .TRUE. ! should DENS_EI_ALL work this way?
        doUncoupledForAddtl = .FALSE.
        call eos_newtonRaphson(vecLen, MODE_DENS_EI, .TRUE., maskPtr,cMask=overallCMask)
        ! Next lines moved from end of eos_newtonRaphson
        eosData(pres+ilo:pres+ihi)=ptotRow(1:rowLen) !Ok
        eosData(gamc+ilo:gamc+ihi)=gamcRow(1:rowLen) !Ok, I guess - if DENS_EI_ALL should work this way? - KW
        eosData(entr+ilo:entr+ihi)=stotRow(1:rowLen) !(Ok, I guess - if DENS_EI_ALL should work this way? - KW)
     end if
     if (mode==MODE_DENS_EI_GATHER .OR. mode==MODE_DENS_EI_RECAL_GATHER) then
        if (ANY(maskPtr((/EOS_PRESELE,EOS_PRESION,EOS_PRESRAD/)) .EQV. .FALSE.)) then
           call Driver_abortFlash("[eos_idealGamma] cannot calculate MODE_DENS_EI_GATHER without component pressure masks.&
                & Set mask appropriately.")
        end if
        eosData(eint+ilo:eint+ihi) = ( &
             eosData(eintIon+ilo:eintIon+ihi)+ &
             eosData(eintEle+ilo:eintEle+ihi)+ &
             eosData(eintRad+ilo:eintRad+ihi))

        if(eos_combinedTempRule==0) then
           call eos_newtonRaphson(vecLen, MODE_DENS_EI, .TRUE., maskPtr,cMask=overallCMask)
           eosData(entr+ilo:entr+ihi)=stotRow(1:rowLen) !(replace with expression based on component states)
        endif
        ! Next lines moved from end of eos_newtonRaphson; pres changed - KW 2010-12-13
        eosData(pres+ilo:pres+ihi) = ( &
             eosData(presIon+ilo:presIon+ihi)+ &
             eosData(presEle+ilo:presEle+ihi)+ &
             eosData(presRad+ilo:presRad+ihi))
#ifdef EOS_FORCE_2T
        eosData(gamc+ilo:gamc+ihi)=eos_gammaEle
#else
        call setCombinedGamc(gamc)
!        eosData(gamc+ilo:gamc+ihi)=gamcRow(1:rowLen) !replace with expression based on component states!
#endif
        call setCombinedTemp()
     end if

    !! Note that we make here a choice for the arbitrary additive constant in the electron entropy 
     if(maskPtr(EOS_ENTRELE)) &
          call setEleEntropy(entrele,tempEle)
     if(maskPtr(EOS_ENTRRAD)) then
        call setRadEntropy(entrrad,tempRad,eintRad)
     end if

!!!! NEW EXPERIMENTAL MODE_DENS_EI_MAT_GATHER
!!!!

  case (MODE_DENS_EI_MAT_GATHER)

     useNRForAddtlCombined = .FALSE.
     doUncoupledForAddtl = .TRUE.
     eosData(presIon+ilo:presIon+ihi) =      eosData(dens+ilo:dens+ihi)  * &
          eosData(eintIon+ilo:eintIon+ihi) * gam1invIon
     eosData(presEle+ilo:presEle+ihi) =      eosData(dens+ilo:dens+ihi)  * &
          eosData(eintEle+ilo:eintEle+ihi) * gam1invEle


     eosData(tempIon+ilo:tempIon+ihi) = eosData(eintIon+ilo:eintIon+ihi) * ggprodinvIon * &
          eosData(abar+ilo:abar+ihi)
     eosData(tempEle+ilo:tempEle+ihi) = eosData(eintEle+ilo:eintEle+ihi) * ggprodinvEle * &
          eosData(abar+ilo:abar+ihi) * Zinv

     eosData(tempRad+ilo:tempRad+ihi) = sqrt(sqrt( eosData(presRad+ilo:presRad+ihi) / eos_ao3 ))


!!     eosData(entr+ilo:entr+ihi) = 0.0

     if (.TRUE.) then
        if (ANY(maskPtr((/EOS_PRESELE,EOS_PRESION,EOS_PRESRAD/)) .EQV. .FALSE.)) then
           call Driver_abortFlash("[eos_idealGamma] cannot calculate MODE_DENS_EI_GATHER without component pressure masks.&
                & Set mask appropriately.")
        end if

        eosData(presRad+ilo:presRad+ihi) = 0.0
        if(eos_combinedTempRule==0) then
           eosData(eint+ilo:eint+ihi) = ( &
                eosData(eintIon+ilo:eintIon+ihi)+ &
                eosData(eintEle+ilo:eintEle+ihi))
           call eos_newtonRaphson(vecLen, MODE_DENS_EI, .TRUE., maskPtr,cMask=(/1,1,0/))
           eosData(entr+ilo:entr+ihi)=stotRow(1:rowLen) !(replace with expression based on component states)
        endif
        eosData(pres+ilo:pres+ihi)       = eosData(presIon+ilo:presIon+ihi)+eosData(presEle+ilo:presEle+ihi)
        call setCombinedGamc(gamc)
        eosData(eint+ilo:eint+ihi) = ( &
             eosData(eintIon+ilo:eintIon+ihi)+ &
             eosData(eintEle+ilo:eintEle+ihi)+ &
             eosData(eintRad+ilo:eintRad+ihi))
        eosData(presRad+ilo:presRad+ihi) = eosData(eintRad+ilo:eintRad+ihi) * eosData(dens+ilo:dens+ihi) / 3.0
        eosData(tempRad+ilo:tempRad+ihi) = sqrt(sqrt( eosData(presRad+ilo:presRad+ihi) / eos_ao3 ))

        eosData(pres+ilo:pres+ihi) = ( &
             eosData(presIon+ilo:presIon+ihi)+ &
             eosData(presEle+ilo:presEle+ihi)+ &
             eosData(presRad+ilo:presRad+ihi))
        call setCombinedTemp()
     end if

    !! Note that we make here a choice for the arbitrary additive constant in the electron entropy 
     if(maskPtr(EOS_ENTRELE)) &
          call setEleEntropy(entrele,tempEle)
     if(maskPtr(EOS_ENTRRAD)) then
        call setRadEntropy(entrrad,tempRad,eintRad)
     end if

!!!! *ADDITIONAL* NEW EXPERIMENTAL MODE_DENS_EI_MAT_EQUI
!!!!

  case (MODE_DENS_EI_MAT_EQUI)

     eosData(tempIon+ilo:tempIon+ihi) = eosData(eintIon+ilo:eintIon+ihi) * ggprodinvIon * &
          eosData(abar+ilo:abar+ihi)
     eosData(tempEle+ilo:tempEle+ihi) = eosData(eintEle+ilo:eintEle+ihi) * ggprodinvEle * &
          eosData(abar+ilo:abar+ihi) * Zinv

     eosData(tempRad+ilo:tempRad+ihi) = sqrt(sqrt( eosData(presRad+ilo:presRad+ihi) / eos_ao3 ))


!!     eosData(entr+ilo:entr+ihi) = 0.0


     useNRForAddtlCombined = .TRUE. ! should DENS_EI_MAT_EQUI work this way?
     doUncoupledForAddtl = .TRUE.   !???????

     if (.TRUE.) then
        if (ANY(maskPtr((/EOS_PRESELE,EOS_PRESION,EOS_PRESRAD/)) .EQV. .FALSE.)) then
           call Driver_abortFlash("[eos_idealGamma] cannot calculate MODE_DENS_EI_GATHER without component pressure masks.&
                & Set mask appropriately.")
        end if
        eosData(presRad+ilo:presRad+ihi) = 0.0
        eosData(eint+ilo:eint+ihi) = ( &
             eosData(eintIon+ilo:eintIon+ihi)+ &
             eosData(eintEle+ilo:eintEle+ihi))

        call eos_newtonRaphson(vecLen, MODE_DENS_EI, .TRUE., maskPtr,cMask=(/1,1,0/))
        eosData(presIon+ilo:presIon+ihi)=pCompRow(EOSCOMP_ION,1:rowLen)
        eosData(presEle+ilo:presEle+ihi)=pCompRow(EOSCOMP_ELE,1:rowLen)
        eosData(eintIon+ilo:eintIon+ihi)=eCompRow(EOSCOMP_ION,1:rowLen)
        eosData(eintEle+ilo:eintEle+ihi)=eCompRow(EOSCOMP_ELE,1:rowLen)
        eosData(gamc+ilo:gamc+ihi)=gamcRow(1:rowLen) ! DENS_EI_MAT_EQUI should work this way. - KW
        eosData(entr+ilo:entr+ihi)=stotRow(1:rowLen) !(replace with expression based on component states)
        eosData(eint+ilo:eint+ihi) = ( &
             eosData(eintIon+ilo:eintIon+ihi)+ &
             eosData(eintEle+ilo:eintEle+ihi)+ &
             eosData(eintRad+ilo:eintRad+ihi))
        eosData(presRad+ilo:presRad+ihi) = eosData(eintRad+ilo:eintRad+ihi) * eosData(dens+ilo:dens+ihi) / 3.0
        eosData(tempRad+ilo:tempRad+ihi) = sqrt(sqrt( eosData(presRad+ilo:presRad+ihi) / eos_ao3 ))

        eosData(tempIon+ilo:tempIon+ihi)=tempRow(1:rowLen)
        eosData(tempEle+ilo:tempEle+ihi)=tempRow(1:rowLen)
        if(eos_combinedTempRule==0) then
        endif

        eosData(pres+ilo:pres+ihi) = ( &
             eosData(presIon+ilo:presIon+ihi)+ &
             eosData(presEle+ilo:presEle+ihi)+ &
             eosData(presRad+ilo:presRad+ihi))
        call setCombinedTemp()
     end if

    !! Note that we make here a choice for the arbitrary additive constant in the electron entropy 
     if(maskPtr(EOS_ENTRELE)) &
          call setEleEntropy(entrele,tempEle)
     if(maskPtr(EOS_ENTRRAD)) then
        call setRadEntropy(entrrad,tempRad,eintRad)
     end if



!!!!
!!!! END OF NEW EXPERIMENTAL MODE_DENS_EI_MAT




!!$  case (MODE_DENS_EI_ALL)
!!$     eosData(pres+ilo:pres+ihi) = eosData(dens+ilo:dens+ihi) * &
!!$          eosData(eint+ilo:eint+ihi) * gam1inv
!!$     eosData(presIon+ilo:presIon+ihi) = (ionFrac * eosData(dens+ilo:dens+ihi)) * &
!!$          eosData(eintIon+ilo:eintIon+ihi) * gam1invIon
!!$     eosData(presEle+ilo:presEle+ihi) = (eleFrac * eosData(dens+ilo:dens+ihi)) * &
!!$          eosData(eintEle+ilo:eintEle+ihi) * gam1invEle
!!$     eosData(presRad+ilo:presRad+ihi) = (radFrac * eosData(dens+ilo:dens+ihi)) * &
!!$          eosData(eintRad+ilo:eintRad+ihi) * gam1invRad
!!$
!!$     eosData(temp+ilo:temp+ihi) = eosData(eint+ilo:eint+ihi) * ggprodinv * &
!!$          eosData(abar+ilo:abar+ihi) * ZpInv
!!$     eosData(tempIon+ilo:tempIon+ihi) = eosData(eintIon+ilo:eintIon+ihi) * ggprodinvIon * &
!!$          eosData(abar+ilo:abar+ihi)
!!$     eosData(tempEle+ilo:tempEle+ihi) = eosData(eintEle+ilo:eintEle+ihi) * ggprodinvEle * &
!!$          eosData(abar+ilo:abar+ihi) * Zinv
!!$     eosData(tempRad+ilo:tempRad+ihi) = eosData(eintRad+ilo:eintRad+ihi) * ggprodinvRad * &
!!$          eosData(abar+ilo:abar+ihi)
!!$
!!$     eosData(entr+ilo:entr+ihi) = 0.0
!!$


  ! density, pressure taken as input
  case (MODE_DENS_PRES)
     ggprodinv = 1. / ggprod
!     gam1inv   = 1. / eos_gammam1 !currently not used
     eosData(eint+ilo:eint+ihi) = eosData(pres+ilo:pres+ihi) * eos_gammam1 / &
                                   eosData(dens+ilo:dens+ihi)
     eosData(temp+ilo:temp+ihi) = eosData(eint+ilo:eint+ihi) * ggprodinv * &
                                   eosData(abar+ilo:abar+ihi) * ZpInv
     eosData(entr+ilo:entr+ihi) = 0.0

!!$     eosData(tempRad+ilo:tempRad+ihi) = sqrt(sqrt( eosData(presRad+ilo:presRad+ihi) / eos_ao3 ))
!!$     eosData(eintRad+ilo:eintRad+ihi)    = 3.0e0 * eosData(presRad+ilo:presRad+ihi) / eosData(dens+ilo:dens+ihi)


  case(MODE_DENS_EI_SELE_GATHER)
     ! We believe eintRad
     eosData(presRad+ilo:presRad+ihi) = eosData(eintRad+ilo:eintRad+ihi) * eosData(dens+ilo:dens+ihi) / 3.0
     eosData(tempRad+ilo:tempRad+ihi) = sqrt(sqrt( eosData(presRad+ilo:presRad+ihi) / eos_ao3 ))
     ! Done with _RAD

     ! We believe entrEle, do not believe eintEle
    select case (seleVariant)
    case(1)
#if(1)
     eosData(tempEle+ilo:tempEle+ihi) = exp( gam1invEle * ( eosData(entrele+ilo:entrele+ihi)/kBoltz &
                                                            + log( eosData(dens+ilo:dens+ihi) ) ) )
#endif
    case(2)
     eosData(tempEle+ilo:tempEle+ihi) = &
          exp( gam1invEle * ( eosData(entrele+ilo:entrele+ihi)/(eos_kergavo*Ye_ARRAYREFERENCE) &
                             + log( eosData(dens+ilo:dens+ihi) * eos_avo**1 *Ye_ARRAYREFERENCE) ) ) &
           * eos_h*eos_hbar / (kBoltz*eos_eMass)                  
    case(3)
     eosData(tempEle+ilo:tempEle+ihi) = &
          exp( gam1invEle * ( eosData(entrele+ilo:entrele+ihi)/(eos_kergavo*Ye_ARRAYREFERENCE) &
                             - (1.0+eos_gammam1Ele) &
                             + log( eosData(dens+ilo:dens+ihi) * eos_avo*Ye_ARRAYREFERENCE) ) ) &
           * eos_h*eos_hbar / (kBoltz*eos_eMass)                  
    case(4)
     eosData(tempEle+ilo:tempEle+ihi) = &
          exp( gam1invEle * ( eosData(entrele+ilo:entrele+ihi)/Ye_ARRAYREFERENCE &
                             + log( eosData(dens+ilo:dens+ihi) * Ye_ARRAYREFERENCE) ) )

    case(5)
     eosData(tempEle+ilo:tempEle+ihi) = exp( gam1invEle * ( eosData(entrele+ilo:entrele+ihi) &
                                                            + log( eosData(dens+ilo:dens+ihi) * Ye_ARRAYREFERENCE) ) )

    case(6)
     eosData(tempEle+ilo:tempEle+ihi) = tempRef / eos_eMassInUAmu * &
                                        exp( gam1invEle * ( eosData(entrele+ilo:entrele+ihi)/Ye_ARRAYREFERENCE &
                                                 + log( eosData(dens+ilo:dens+ihi)*densRefInv  * Ye_ARRAYREFERENCE) ) )
    case(7)
#if(1)
     eosData(tempEle+ilo:tempEle+ihi) = (exp( eosData(entrele+ilo:entrele+ihi)) &
                                                            * eosData(dens+ilo:dens+ihi) )**gam1invEle
#endif
    case(8)
     eosData(tempEle+ilo:tempEle+ihi) = (     eosData(entrele+ilo:entrele+ihi)  &
                                                            * eosData(dens+ilo:dens+ihi) )**gam1invEle
     end select
     eosData(eintEle+ilo:eintEle+ihi) = ggprodEle * eosData(tempEle+ilo:tempEle+ihi) &
          / eosData(abar+ilo:abar+ihi) * dynamicZ
     forall (i=ilo:ihi, (eosData(eintEle+i)<eos_smallEele))
        eosData(eintEle+i) = eos_smallEele
        eosData(tempEle+i) = eos_smallEele  * ggprodinvEle * eosData(abar+i) * Zinv
     end forall
     eosData(presEle+ilo:presEle+ihi) = eos_gasConstant * &
          eosData(dens+ilo:dens+ihi) * &
          eosData(tempEle+ilo:tempEle+ihi) / eosData(abar+ilo:abar+ihi) * dynamicZ
     ! eintEle and other _ELE have now been recomputed from entrEle

     if (ANY(maskPtr((/EOS_PRESELE,EOS_PRESION,EOS_PRESRAD/)) .EQV. .FALSE.)) then
        call Driver_abortFlash("[eos_idealGamma] cannot calculate MODE_DENS_EI_SELE_GATHER without component pressure masks.&
             & Set mask appropriately.")
     end if

     ! We believe eint (combined), do NOT believe eintIon
     eosData(eintIon+ilo:eintIon+ihi) = max( eos_smallEion, &
                                            eosData(eint+ilo:eint+ihi) - &
                                            eosData(eintEle+ilo:eintEle+ihi) - &
                                            eosData(eintRad+ilo:eintRad+ihi)  )
     ! eintIon has now been recomputed

     eosData(presIon+ilo:presIon+ihi) = (1 * eosData(dens+ilo:dens+ihi)) * &
          eosData(eintIon+ilo:eintIon+ihi) * gam1invIon
     eosData(tempIon+ilo:tempIon+ihi) = eosData(eintIon+ilo:eintIon+ihi) * ggprodinvIon * &
          eosData(abar+ilo:abar+ihi)
     ! Other _ION have now been recomputed from eintIon

#if(0)
     eosData(eint+ilo:eint+ihi) = ( &
          eosData(eintIon+ilo:eintIon+ihi)+ &
          eosData(eintEle+ilo:eintEle+ihi)+ &
          eosData(eintRad+ilo:eintRad+ihi))
#endif
     useNRForAddtlCombined = .FALSE.
     doUncoupledForAddtl = .TRUE.

     if(eos_combinedTempRule==0) then
     ! Do Newton-Raphson iteration, which calls eos_byTemp repeatedly.
     ! In this case, we call it mostly for the effect of computing an
     ! effective combined temperature for TEMP_VAR.  PRES_VAR is also
     ! set by this call.
        call eos_newtonRaphson(vecLen, MODE_DENS_EI, .TRUE., maskPtr,cMask=overallCMask)
        eosData(entr+ilo:entr+ihi)=stotRow(1:rowLen) !(replace with expression based on component states)
     end if
     
     ! Next lines moved from end of eos_newtonRaphson; pres changed - KW 2010-12-13
     eosData(pres+ilo:pres+ihi) = ( &
          eosData(presIon+ilo:presIon+ihi)+ &
          eosData(presEle+ilo:presEle+ihi)+ &
          eosData(presRad+ilo:presRad+ihi))
#ifdef EOS_FORCE_2T
     eosData(gamc+ilo:gamc+ihi)=eos_gammaEle
#else
     call setCombinedGamc(gamc)
!     eosData(gamc+ilo:gamc+ihi)=gamcRow(1:rowLen) !replace with expression based on component states!
#endif
     call setCombinedTemp()

  ! unrecognized value for mode
  case default
     call Driver_abortFlash("[eos_idealGamma] Unrecognized input mode given to Eos")
  end select



  if(present(mask)) then

     if(mask(EOS_DPT)) then
        dpt = (EOS_DPT-1)*vecLen
        if(useNRForAddtlCombined) then
           eosData(dpt+ilo:dpt+ihi) = dptRow(1:rowLen)
        else if (doUncoupledForAddtl) then
           eosData(dpt+ilo:dpt+ihi) = (  eos_gasConstant * eosData(dens+ilo:dens+ihi) * &
                                        ( eosData(tempIon+ilo:tempIon+ihi) + eosData(tempEle+ilo:tempEle+ihi)*dynamicZ)&
                                         / eosData(abar+ilo:abar+ihi)  + &
                                        4*eos_ao3 ) / eosData(tempEle+ilo:tempEle+ihi)
        else
           eosData(dpt+ilo:dpt+ihi) = eos_gasConstant*eosData(dens+ilo:dens+ihi) / eosData(abar+ilo:abar+ihi)
        end if
     end if
     if(mask(EOS_DPD)) then
        dpd = (EOS_DPD-1)*vecLen
        eosData(dpd+ilo:dpd+ihi) = eos_gasConstant*eosData(temp+ilo:temp+ihi) / eosData(abar+ilo:abar+ihi)
     end if
     if(mask(EOS_DET))then
        det = (EOS_DET-1)*vecLen
        if(useNRForAddtlCombined) then
           eosData(det+ilo:det+ihi) = detRow(1:rowLen)
        else if (doUncoupledForAddtl) then
#ifdef EOS_FORCE_2T
           call averageVectRByPressureFraction(det, &
                eosData(eintIon+ilo:eintIon+ihi) / eosData(tempIon+ilo:tempIon+ihi), &
                eosData(eintEle+ilo:eintEle+ihi) / eosData(tempEle+ilo:tempEle+ihi))
#else
           call averageVectRByPressureFraction(det, &
                eosData(eintIon+ilo:eintIon+ihi) / eosData(tempIon+ilo:tempIon+ihi), &
                eosData(eintEle+ilo:eintEle+ihi) / eosData(tempEle+ilo:tempEle+ihi), &
                4 * eosData(eintRad+ilo:eintRad+ihi) / eosData(tempRad+ilo:tempRad+ihi))
#endif
        else
           eosData(det+ilo:det+ihi) = ggprod * Zp / eosData(abar+ilo:abar+ihi)
        end if
     end if
     if(mask(EOS_DED))then 
        ded = (EOS_DED-1)*vecLen
        eosData(ded+ilo:ded+ihi) = 0.
     end if

    ! Entropy derivatives   
     if (mask(EOS_DST)) then
        if (mask(EOS_DET) .AND. mask(EOS_DPT)) then
           det = (EOS_DET-1)*vecLen
           dpt = (EOS_DPT-1)*vecLen
           dst = (EOS_DST-1)*vecLen
           eosData(dst+ilo:dst+ihi) = ( (eosData(dpt+ilo:dpt+ihi)  / eosData(dens+ilo:dens+ihi) + eosData(det+ilo:det+ihi)) -&
                &                      (eosData(pres+ilo:pres+ihi)/ eosData(dens+ilo:dens+ihi) + eosData(eint+ilo:eint+ihi))/ &
                &                      eosData(temp+ilo:temp+ihi) ) / eosData(temp+ilo:temp+ihi)
        else
           call Driver_abortFlash("[eos_idealGamma] Cannot calculate EOS_DST without EOS_DET and EOS_DPT")
        end if
     end if
     if (mask(EOS_DSD)) then
        if (mask(EOS_DED) .AND. mask(EOS_DPD)) then
           dsd = (EOS_DSD-1)*vecLen
           ded = (EOS_DED-1)*vecLen
           dpd = (EOS_DPD-1)*vecLen
           eosData(dsd+ilo:dsd+ihi) = &
               ( ((eosData(dpd+ilo:dpd+ihi) - eosData(pres+ilo:pres+ihi)/eosData(dens+ilo:dens+ihi)) / &
        &          eosData(dens+ilo:dens+ihi)) + eosData(ded+ilo:ded+ihi)) / eosData(temp+ilo:temp+ihi)
        else
           call Driver_abortFlash("[eos_idealGamma] Cannot calculate EOS_DSD without EOS_DED and EOS_DPD")
        end if
     end if


     if(mask(EOS_PEL))then 
        pel = (EOS_PEL-1)*vecLen
        eosData(pel+ilo:pel+ihi) = 0.
     end if
     if(mask(EOS_NE))then 
        ne = (EOS_NE-1)*vecLen
        eosData(ne+ilo:ne+ihi) = 0.
     end if
     if(mask(EOS_ETA))then 
        call Driver_abortFlash("[eos_idealGamma] cannot calculate ETA in the multiTemp / Gamma implementation.")
        eta = (EOS_ETA-1)*vecLen
        eosData(eta+ilo:eta+ihi) = 0.
     end if
     
     if(mask(EOS_CV))then
        if(mask(EOS_DET)) then
           c_v = (EOS_CV-1)*vecLen
           if(useNRForAddtlCombined) then
              eosData(c_v+ilo:c_v+ihi) = cvRow(1:rowLen)
           else if (doUncoupledForAddtl) then
              eosData(c_v+ilo:c_v+ihi) = eosData(det+ilo:det+ihi)
           else
              eosData(c_v+ilo:c_v+ihi) = eosData(det+ilo:det+ihi)
           end if
        else
           call Driver_abortFlash("[eos_idealGamma] cannot calculate C_V without DET.  Set mask appropriately.")
        end if
     end if
     ! ideal gas -- all gammas are equal
     if(mask(EOS_CP))then
        if(mask(EOS_CV).and.mask(EOS_DET)) then
           c_p = (EOS_CP-1)*vecLen
           if(useNRForAddtlCombined) then
              eosData(c_p+ilo:c_p+ihi) = cpRow(1:rowLen)
           else if (doUncoupledForAddtl) then
              call averageScalByPressureFraction(c_p, &
                   eos_gammaIon, eos_gammaEle, eos_gammaRad)
              eosData(c_p+ilo:c_p+ihi) = eosData(c_p+ilo:c_p+ihi)*eosData(c_v+ilo:c_v+ihi)
           else
              eosData(c_p+ilo:c_p+ihi) = eos_gamma*eosData(c_v+ilo:c_v+ihi)
           end if
        else
           call Driver_abortFlash("[eos_idealGamma] cannot calculate C_P without C_V and DET.  Set mask appropriately.")
        end if
        
     end if

     if(mask(EOS_CVION))then
           c_v = (EOS_CVION-1)*vecLen
           eosData(c_v+ilo:c_v+ihi) = ggprod  / eosData(abar+ilo:abar+ihi)
     end if
     if(mask(EOS_CVELE))then
           c_v = (EOS_CVELE-1)*vecLen
           eosData(c_v+ilo:c_v+ihi) = ggprod * dynamicZ / eosData(abar+ilo:abar+ihi)
     end if
  end if

#ifdef DEBUG_MAR2012
  if (ANY(abs(eosData(eint+ilo:eint+ihi) - ( &
       eosData(eintIon+ilo:eintIon+ihi)+ &
       eosData(eintEle+ilo:eintEle+ihi)+ &
       eosData(eintRad+ilo:eintRad+ihi))) .GE. 1.e-10*abs(eosData(eint+ilo:eint+ihi)))) then
     print*,'EINT2',eosData(eint+ilo:eint+ihi)
     print*,'ESUM2',(eosData(eintIon+ilo:eintIon+ihi)+eosData(eintEle+ilo:eintEle+ihi))
     print*,'EION2',eosData(eintIon+ilo:eintIon+ihi)
     print*,'EELE2',eosData(eintEle+ilo:eintEle+ihi)
  end if
#endif

  if (arrayBoundHackMode .eqv. .true.) then
     deallocate(maskPtr)
  end if


  !! Close up arrays if previously allocated
#ifndef FIXEDBLOCKSIZE  
  call eos_vecDealloc()
#endif

contains
  subroutine setEleEntropy(entrComponent,tempComponent)
    integer,intent(IN) :: entrComponent, tempComponent
    real :: prefactor

    select case (seleVariant)
    case(1)
       prefactor = kBoltz
    case(2,3)
       prefactor = eos_kergavo * Ye
    case(4,6)
       prefactor = Ye
    case(5,7)
       prefactor = 1.0
    case default
       ! call Driver_evolveFlash('eos_multiTypeByTemp - unsupported eos_entrEleScaleChoice value')
    end select
    select case (seleVariant)
    case(1)
    eosData(entrComponent+ilo:entrComponent+ihi) = &
         kBoltz * &
         (eos_gammam1Ele * log( eosData(tempComponent+ilo:tempComponent+ihi) ) &
                      -    log( eosData(dens+ilo:dens+ihi)) )
    case(2)
    eosData(entrComponent+ilo:entrComponent+ihi) = &
         eos_avo     * eosData(zbar+ilo:zbar+ihi)/eosData(abar+ilo:abar+ihi) * kBoltz * &
         (eos_gammam1Ele * log( eosData(tempComponent+ilo:tempComponent+ihi) * &
                                kBoltz *eos_eMass / (eos_h*eos_hbar) ) &
                      -    log( eosData(dens+ilo:dens+ihi) * eos_avo**1  *eosData(zbar+ilo:zbar+ihi) / eosData(abar+ilo:abar+ihi)) )
    case(3)
    eosData(entrComponent+ilo:entrComponent+ihi) = &
         eos_kergavo * Ye_ARRAYREFERENCE * &
         (1.0 +  &
          eos_gammam1Ele * (1.0 + log( eosData(tempComponent+ilo:tempComponent+ihi) * &
                                       kBoltz *eos_eMass / (eos_h*eos_hbar) ) ) &
                      -    log( eosData(dens+ilo:dens+ihi) * eos_avo * Ye_ARRAYREFERENCE) )
    case(4)
    eosData(entrComponent+ilo:entrComponent+ihi) = &
         Ye_ARRAYREFERENCE * &
         (eos_gammam1Ele * log( eosData(tempComponent+ilo:tempComponent+ihi)) &
                      -    log( eosData(dens+ilo:dens+ihi) * Ye_ARRAYREFERENCE) )
    case(5)
    eosData(entrComponent+ilo:entrComponent+ihi) = &
         (eos_gammam1Ele * log( eosData(tempComponent+ilo:tempComponent+ihi)) &
                      -    log( eosData(dens+ilo:dens+ihi) * Ye_ARRAYREFERENCE) )
    case(6)
    eosData(entrComponent+ilo:entrComponent+ihi) = &
         Ye_ARRAYREFERENCE * &
         log( (eosData(tempComponent+ilo:tempComponent+ihi)*tempRefInv*eos_eMassInUAmu)**eos_gammam1Ele &
                  / ( eosData(dens+ilo:dens+ihi)*densRefInv  * Ye_ARRAYREFERENCE) )
    case(7)
    eosData(entrComponent+ilo:entrComponent+ihi) = &
         (                 log( eosData(tempComponent+ilo:tempComponent+ihi)**eos_gammam1Ele &
                           /    eosData(dens+ilo:dens+ihi)                    ) )
    case(8)
    eosData(entrComponent+ilo:entrComponent+ihi) = &
                                eosData(tempComponent+ilo:tempComponent+ihi)**eos_gammam1Ele &
                           /    eosData(dens+ilo:dens+ihi)                    
    end select
  end subroutine setEleEntropy

  subroutine setRadEntropy(entrComponent,tempComponent,eintComponent)
    integer,intent(IN) :: entrComponent, tempComponent, eintComponent

#define sradVariant 4
    select case (sradVariant)
    case(1)
    eosData(entrComponent+ilo:entrComponent+ihi) = &
         kBoltz * &
         (eos_gammam1Rad * log( eosData(tempComponent+ilo:tempComponent+ihi) ) &
                      -    log( eosData(dens+ilo:dens+ihi)) )
    case(2)
    eosData(entrComponent+ilo:entrComponent+ihi) = &
         eos_kergavo * eosData(zbar+ilo:zbar+ihi)/eosData(abar+ilo:abar+ihi) * kBoltz * &
         (eos_gammam1Rad * log( eosData(tempComponent+ilo:tempComponent+ihi) * &
                                kBoltz *eos_eMass / (eos_h*eos_hbar) ) &
                      -    log( eosData(dens+ilo:dens+ihi) * eos_avo *eosData(zbar+ilo:zbar+ihi) / eosData(abar+ilo:abar+ihi)) )
    case(3)
    eosData(entrComponent+ilo:entrComponent+ihi) = &
         eos_kergavo * Ye_ARRAYREFERENCE * kBoltz * &
         (eos_gammam1Rad * log( eosData(tempComponent+ilo:tempComponent+ihi) * &
                                kBoltz *eos_eMass / (eos_h*eos_hbar) ) &
                      -    log( eosData(dens+ilo:dens+ihi) * eos_avo * Ye_ARRAYREFERENCE) )
    case(4)
       where ( eosData(tempComponent+ilo:tempComponent+ihi) > 0.0 )
          eosData(entrComponent+ilo:entrComponent+ihi) = &
               eos_gammaRad * eosData(eintComponent+ilo:eintComponent+ihi) / eosData(tempComponent+ilo:tempComponent+ihi)
       elsewhere
          eosData(entrComponent+ilo:entrComponent+ihi) = 0.0
       end where
    case(5)
    eosData(entrComponent+ilo:entrComponent+ihi) = &
         (eos_gammam1Rad * log( eosData(tempComponent+ilo:tempComponent+ihi)) &
                      -    log( eosData(dens+ilo:dens+ihi) * Ye_ARRAYREFERENCE) )
    case(6)
    eosData(entrComponent+ilo:entrComponent+ihi) = &
         (                 log( eosData(tempComponent+ilo:tempComponent+ihi)**eos_gammam1Rad &
                           /    eosData(dens+ilo:dens+ihi)                    ) )
    case(7)
    eosData(entrComponent+ilo:entrComponent+ihi) = &
         (                 log( eosData(tempComponent+ilo:tempComponent+ihi)**eos_gammam1Rad &
                           /    eosData(dens+ilo:dens+ihi)                    ) )
    case(8)
    eosData(entrComponent+ilo:entrComponent+ihi) = &
                                eosData(tempComponent+ilo:tempComponent+ihi)**eos_gammam1Rad &
                           /    eosData(dens+ilo:dens+ihi)                    
    end select
  end subroutine setRadEntropy

  subroutine setCombinedGamc(gamcCombined)
    integer, intent(in) :: gamcCombined

    eosData(gamcCombined+ilo:gamcCombined+ihi) = &
         ( eos_gammaIon*eosData(presIon+ilo:presIon+ihi)+ &
           eos_gammaEle*eosData(presEle+ilo:presEle+ihi)+ &
           eos_gammaRad*eosData(presRad+ilo:presRad+ihi) ) / eosData(pres+ilo:pres+ihi)

  end subroutine setCombinedGamc

  subroutine setCombinedTemp()
    select case (eos_combinedTempRule)
    case(1)
       eosData(temp+ilo:temp+ihi) = eosData(tempIon+ilo:tempIon+ihi)
    case(2)
       eosData(temp+ilo:temp+ihi) = eosData(tempEle+ilo:tempEle+ihi)
    case(3)
       eosData(temp+ilo:temp+ihi) = eosData(tempRad+ilo:tempRad+ihi)

       ! otherwise: either already done, or do not care.
    end select

  end subroutine setCombinedTemp

  subroutine averageScalByPressureFraction(combined,scal1,scal2,scal3)
    integer, intent(in) :: combined
    real, intent(in) :: scal1,scal2,scal3

    eosData(combined+ilo:combined+ihi) = &
         ( scal1*eosData(presIon+ilo:presIon+ihi)+ &
           scal2*eosData(presEle+ilo:presEle+ihi)+ &
           scal3*eosData(presRad+ilo:presRad+ihi) ) / eosData(pres+ilo:pres+ihi)

  end subroutine averageScalByPressureFraction

  subroutine averageVectByPressureFraction(combined,compIon,compEle,compRad)
    integer, intent(in) :: combined,compIon,compEle,compRad

    eosData(combined+ilo:combined+ihi) = &
         ( eosData(compIon+ilo:compIon+ihi)*eosData(presIon+ilo:presIon+ihi)+ &
           eosData(compEle+ilo:compEle+ihi)*eosData(presEle+ilo:presEle+ihi)+ &
           eosData(compRad+ilo:compRad+ihi)*eosData(presRad+ilo:presRad+ihi) ) / eosData(pres+ilo:pres+ihi)

  end subroutine averageVectByPressureFraction

  subroutine averageVectRByPressureFraction(combined,compIon,compEle,compRad)
    integer, intent(in) :: combined
    real, intent(in),dimension(1:rowLen) :: compIon,compEle
    real, intent(in),dimension(1:rowLen),OPTIONAL :: compRad

#ifdef EOS_FORCE_2T
    eosData(combined+ilo:combined+ihi) = &
         ( compIon(1:rowLen)*eosData(presIon+ilo:presIon+ihi)+ &
           compEle(1:rowLen)*eosData(presEle+ilo:presEle+ihi) ) / eosData(pres+ilo:pres+ihi)
#else
    eosData(combined+ilo:combined+ihi) = &
         ( compIon(1:rowLen)*eosData(presIon+ilo:presIon+ihi)+ &
           compEle(1:rowLen)*eosData(presEle+ilo:presEle+ihi)+ &
           compRad(1:rowLen)*eosData(presRad+ilo:presRad+ihi) ) / eosData(pres+ilo:pres+ihi)
#endif

  end subroutine averageVectRByPressureFraction


  subroutine eos_newtonRaphson(vecLen, mode, enerWanted, eosMask,cMask)
  use eos_vecData, ONLY:  denRow, etotRow, abarRow, zbarRow, &
       deaRow, dezRow, dsdRow, dstRow, &
       tempRadRow,tempIonRow,tempEleRow
    implicit none

    integer,intent(IN) :: vecLen, mode
    logical,intent(IN) :: enerWanted

    logical,optional, dimension(EOS_VARS+1:EOS_NUM),INTENT(in)::eosMask
    integer,optional, dimension(EOSCOMP_NUM_COMPONENTS),INTENT(in)::cMask

    integer :: k, kk
    integer :: rowBegin,rowEnd
  ! declare some local storage for the results of the Newton iteration
    real,dimension(rowLen)::  ewantRow, tnew, error,pwantRow
  !  local storage for forcedConstantInput -- could be allocatable, but might be so slow
  !  that it's not worth the small storage save.
    real,dimension(rowLen)::  psaveRow, esaveRow

    rowBegin = 1
    rowEnd = rowLen

    do k = 1,rowLen
       kk = k+ilo-1
       tempRow(k)    = eosData(temp+kk)
       denRow(k)     = eosData(dens+kk)
       ewantRow(k)   = eosData(eint+kk)   ! store desired internal energy for MODE_DENS_EI etc. cases
       pwantRow(k)   = eosData(pres+kk)   ! store desired pressure for MODE_DENS_PRES etc. cases
       abarRow(k)    = eosData(abar+kk)   ! NOTE user is assumed to know what they're doing!
       zbarRow(k)    = eosData(zbar+kk)   ! NOTE compare to normal Helmholtz where we don't assume
    end do

    !! Save the input parameters if eos_forceConstantInput is defined
    if (eos_forceConstantInput) then
       esaveRow = ewantRow
       psaveRow = pwantRow
    end if

    ! Initialize the errors
    error(:) = 0.0e0

    ! Do the first eos call with all the zones in the pipe
    !  NOTE that eos_byTemp can ONLY operate in the equivalent of
    !  MODE_DENS_TEMP, as it returns pressure, energy and derivatives only
    !  So if you send in a crappy temperature here, you'll get a crappy starting
    !  position and the iteration won't converge.
    !  Initial temperature here is what is stored in the grid, even though we 
    !    SUSPECT this is not in equilibrium (or we wouldn't be calling Eos if it was fine)
    if (.NOT. all(tempRow(1:rowLen) > 0)) &
         print*,'eos_nr WARN:',tempRow(1:rowLen)

    tempRow(1:rowLen) = max(0.0,eos_smallT,tempRow(1:rowLen))
    if (.NOT. present(cMask)) then
       tempRadRow(1:rowLen) = tempRow(1:rowLen)
       tempIonRow(1:rowLen) = tempRow(1:rowLen)
       tempEleRow(1:rowLen) = tempRow(1:rowLen)
    else
       if (cMask(EOSCOMP_RAD).NE.0) tempRadRow(1:rowLen) = tempRow(1:rowLen)
       if (cMask(EOSCOMP_ION)+cMask(EOSCOMP_ELE).NE.0) tempIonRow(1:rowLen) = tempRow(1:rowLen)
       if (cMask(EOSCOMP_ELE).NE.0) tempEleRow(1:rowLen) = tempRow(1:rowLen)
    end if
    call eos_byTemp(rowBegin,rowEnd,eosMask,componentMask=cMask,eleFrac=eleFrac,ggprodEle=ggprodEle)
    !  Now eos_byTemp has returned ptotRow, etotRow, dXXtRow, and gamcRow


    !  Create initial condition
    do k = rowBegin, rowEnd

       if (enerWanted) then
          !  ewantRow is our desired EI input
          tnew(k) = tempRow(k) - (etotRow(k) - ewantRow(k)) / detRow(k)
       else                    !pressure is wanted
          tnew(k) = tempRow(k) - (ptotRow(k) - pwantRow(k)) / dptRow(k)
       end if

       ! Don't allow the temperature to change by more than an order of magnitude 
       ! in a single iteration
       if (tnew(k) .GT. 10.e0*tempRow(k)) tnew(k) =  & 
            &           10.e0*tempRow(k)
       if (tnew(k) .LT. 0.1e0*tempRow(k)) tnew(k) =  & 
            &           0.1e0*tempRow(k)

       ! Compute the error
       error(k) = abs((tnew(k) - tempRow(k)) / tempRow(k))

       ! Store the new temperature
       tempRow(k) = tnew(k)

       ! Check if we are freezing, if so set the temperature to smallt, and adjust 
       ! the error so we don't wait for this one
       if (tempRow(k) .LT. eos_smallt) then
          tempRow(k) = eos_smallt
          error(k)    = 0.1*eos_tol
       endif

    enddo

    ! Loop over the zones individually now
    do k = rowBegin, rowEnd
       do i = 2, eos_maxNewton
          if (error(k) .LT. eos_tol) goto 70

          if (.NOT. present(cMask)) then
             tempRadRow(k) = tempRow(k)
             tempIonRow(k) = tempRow(k)
             tempEleRow(k) = tempRow(k)
          else
             if (cMask(EOSCOMP_RAD).NE.0) tempRadRow(k) = tempRow(k)
             if (cMask(EOSCOMP_ION).NE.0) tempIonRow(k) = tempRow(k)
             if (cMask(EOSCOMP_ELE).NE.0) tempEleRow(k) = tempRow(k)
          end if
          ! do eos only over this single item
          call eos_byTemp(k,k,eosMask,componentMask=cMask,eleFrac=eleFrac,ggprodEle=ggprodEle)

          if (enerWanted) then
             tnew(k) = tempRow(k) - (etotRow(k) - ewantRow(k)) / detRow(k)
          else                    !pressure is wanted
             tnew(k) = tempRow(k) - (ptotRow(k) - pwantRow(k)) / dptRow(k)
          end if

          ! Don't allow the temperature to change by more than an order of magnitude 
          ! in a single iteration
          if (tnew(k) .GT. 10.e0*tempRow(k)) tnew(k) =  & 
               &              10.e0*tempRow(k)
          if (tnew(k) .LT. 0.1e0*tempRow(k)) tnew(k) =  & 
               &              0.1e0*tempRow(k)

          ! Compute the error
          error(k) = abs((tnew(k) - tempRow(k)) / tempRow(k))

          ! Store the new temperature
          tempRow(k) = tnew(k)

          ! Check if we are freezing, if so set the temperature to eos_smallt, and adjust 
          ! the error so we don't wait for this one
          if (tempRow(k) .LT. eos_smallt) then
             tempRow(k) = eos_smallt
             error(k)    = .1*eos_tol
          endif

       end do  ! end of Newton iterations loop.  Failure drops below, success goes to 70

       ! Land here if too many iterations are needed -- failure

       print *, ' '
       print *, 'Newton-Raphson failed in routine multiTemp idealGamma Eos'

       select case(mode)
       case(MODE_DENS_EI)
          print *, '(e and rho as input):'
       case(MODE_DENS_PRES)
          print *, '(pres and rho as input):'
       case default
          print *, 'Unusual Eos mode:', mode
       end select
       print *, ' '
       print *, 'too many iterations'
       print *, ' '
       print *, ' k    = ', k,rowBegin,rowEnd,vecLen
       print *, ' temp = ', tempRow(k)
       print *, ' dens = ', denRow(k)
       print *, ' pres = ', ptotRow(k)

       call Driver_abortFlash('[eos_idealGamma] Error: too many iteration in Newton-Raphson')


       ! Land here if the Newton iteration converged
       !  jumps out of the iterations, but then continues to the next vector location

70     continue           
    end do

    ! Crank through the entire eos one last time

    if (.NOT. present(cMask)) then
       tempRadRow(1:rowLen) = tempRow(1:rowLen)
       tempIonRow(1:rowLen) = tempRow(1:rowLen)
       tempEleRow(1:rowLen) = tempRow(1:rowLen)
    else
       if (cMask(EOSCOMP_RAD).NE.0) tempRadRow(1:rowLen) = tempRow(1:rowLen)
       if (cMask(EOSCOMP_ION).NE.0) tempIonRow(1:rowLen) = tempRow(1:rowLen)
       if (cMask(EOSCOMP_ELE).NE.0) tempEleRow(1:rowLen) = tempRow(1:rowLen)
    end if
    call eos_byTemp(rowBegin,rowEnd,eosMask,componentMask=cMask,eleFrac=eleFrac,ggprodEle=ggprodEle)

    ! Fill the FLASH arrays with the results.  
    !  In MODE_DENS_EI, we should be generating temperature and pressure (plus gamma and entropy)

    select case(mode)
    case(MODE_DENS_EI)
       eosData(temp+ilo:temp+ihi)=tempRow(1:rowLen)
!       eosData(pres+ilo:pres+ihi)=ptotRow(1:rowLen) ! moved up to callers - KW 2010-12-13
       !  Update the energy to be the true energy, instead of the energy we were trying to meet
       !  ConstantInput LBR and KW believe this is WRONG -- the input arrays should not be changed
       if (eos_forceConstantInput)  then
          eosData(eint+ilo:eint+ihi) = esaveRow(1:rowLen)
       else
          eosData(eint+ilo:eint+ihi) = etotRow(1:rowLen)
       end if
    case(MODE_DENS_PRES)
       eosData(temp+ilo:temp+ihi)=tempRow(1:rowLen)
       ! Update the pressure to be the equilibrium pressure, instead of the pressure we were trying to meet
       !  ConstantInput LBR and KW believe this is wrong.  See notes at the top of the routine
       if (eos_forceConstantInput) then
          eosData(pres+ilo:pres+ihi) = psaveRow(1:rowLen)
       else
          eosData(pres+ilo:pres+ihi) = ptotRow(1:rowLen)
       end if
       eosData(eint+ilo:eint+ihi)=etotRow(1:rowLen)
    end select

!    eosData(gamc+ilo:gamc+ihi)=gamcRow(1:rowLen) ! moved up to callers - KW 2010-12-13
!    eosData(entr+ilo:entr+ihi)=stotRow(1:rowLen) ! moved up to callers - KW 2010-12-13
#ifdef DEBUG_DEC2010
    print*,'NR+I',eCompRow(EOSCOMP_ION,1:rowLen)
!    print*,'NR+e',eCompRow(EOSCOMP_ELE,1:rowLen)
!    print*,'NR+R',eCompRow(EOSCOMP_RAD,1:rowLen)
#endif

  end subroutine eos_newtonRaphson


end subroutine eos_idealGamma

!! FOR FUTURE  : This section is not in use in FLASH 3 yet. none
!! of the current setups use entropy. This will be taken care of 
!! in future releases

!!..no matter what the input mode compute the entropy
!!..ignore the -chemical_potential*number_density part for now
!!$  dens_inv = 1.0e0/eosData(dens+1:+vecLen)
!!$  temp_inv = 1.0e0/eosData(temp+1:+vecLen)
!!$  stot     = (pres*dens_inv + eosData(eint+1:+vecLen))*temp_inv 
!!$  dstotdd  = (eosData(EOS_DPD)*dens_inv - pres*dens_inv*dens_inv + eosData(EOS_DED))*temp_inv
!!$  dstotdt  = (eosData(EOS_DPT)*dens_inv + eosData(EOS_DET))*temp_inv  - (pres*dens_inv + eosData(eint+1:+vecLen)) * temp_inv*temp_inv 
!!$  

