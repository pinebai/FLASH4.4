!!****if* source/physics/Eos/EosMain/multiTemp/Multitype/Eos
!!
!!
!! NAME
!!
!!  Eos
!!
!! SYNOPSIS
!!
!!       call Eos(integer(IN) :: mode,
!!                integer(IN) :: vecLen,
!!                real(INOUT) :: eosData(vecLen*EOS_NUM),
!!      optional, real(IN)    :: massFrac(vecLen*NSPECIES),
!!      optional, logical(IN),target :: mask(EOS_VARS+1:EOS_NUM),
!!      optional, integer(IN) :: vecBegin,
!!      optional, integer(IN) :: vecEnd    )
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
!!  vecBegin : Index of first cell in eosData to handle.
!!             Can be used to limit operation to a subrange of cells, untested.
!!             If not present, the default is 1.
!!  vecEnd   : Index of last cell in eosData to handle.
!!             Can be used to limit operation to a subrange of cells, untested.
!!             If not present, the default is vecLen.
!!
!! PARAMETERS
!!
!!     gamma        :  Ratio of specific heats for the simulated gas
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

#ifdef DEBUG_ALL
#define DEBUG_EOS
#endif
#define ORIGINAL_GAMC_AVERAGE

#define WANT_ENER 2
#define WANT_PRES 3
#define WANT_ENTR 4

subroutine Eos(mode, vecLen, eosData, massFrac, mask, vecBegin,vecEnd, diagFlag)

!==============================================================================
  use Eos_data, ONLY : eos_gasConstant, &
       eos_smallT, eos_tol, eos_maxNewton, &
       eos_meshMe, &
       eos_combinedTempRule
  use eos_helmConstData, ONLY: eos_ao3, eos_kerg
  use eos_mgammaData, ONLY: eos_gammam1j,  eos_ggprodj, eos_gc, &
       eos_gammaEle, &
       eos_gammam1Ele, eos_gammam1Rad, &
       eos_maxFactorUp, eos_maxFactorDown
  use Logfile_interface, ONLY : Logfile_stamp
  use Driver_interface, ONLY : Driver_abortFlash
  use Multispecies_interface, ONLY: Multispecies_getSumInv, Multispecies_getSumFrac
  use Hydro_interface, ONLY : Hydro_recalibrateEints

  use eos_mtInterface, ONLY : eos_multiTypeByTemp

  use eos_vecData, ONLY:  stotRow, eCompRow, &
       dpdRow, dedRow, pelRow, neRow, etaRow, cvRow, cpRow

  use Simulation_interface, ONLY: Simulation_mapIntToStr
#include "Eos_components.h"
  implicit none

#include "constants.h"
#include "Eos.h"
#include "Flash.h"
#include "Multispecies.h"

  !     Arguments
  integer, INTENT(in) :: mode, vecLen
  real,INTENT(inout), dimension(EOS_NUM*vecLen) :: eosData 
  real, optional, INTENT(in),dimension(NSPECIES*vecLen)    :: massFrac
  logical, optional, INTENT(in),target,dimension(EOS_VARS+1:EOS_NUM) :: mask
  integer,optional,INTENT(in) :: vecBegin,vecEnd
  integer, optional, INTENT(out)    :: diagFlag

  interface
     subroutine eos_newtonRaphson(vecLen, eosData, ilo,ihi, mode, whatWanted, &
          massFrac, eosMaskPlus,cMask,forceConstantInput)
       implicit none
       integer,intent(IN) :: vecLen, mode, ilo,ihi
       real,INTENT(inout), dimension(EOS_NUM*vecLen) :: eosData 
       integer,intent(IN) :: whatWanted
       real,   INTENT(in),dimension(NSPECIES*vecLen) :: massFrac
       logical,optional, dimension(EOS_VARS+1:EOS_NUM),INTENT(INOUT)::eosMaskPlus
       integer,optional, dimension(EOSCOMP_NUM_COMPONENTS),INTENT(in)::cMask
       logical,optional,intent(in) :: forceConstantInput
     end subroutine eos_newtonRaphson
  end interface

  ! This is the variable that is used internally -- set to false unless mask comes in
  logical,save,target, dimension(EOS_VARS+1:EOS_NUM) :: maskInternal
  logical,dimension(EOS_VARS+1:EOS_NUM) :: maskPlus
  data maskInternal / EOS_DERIVS * .FALSE. /
  logical,pointer, dimension(:) :: maskPtr

#ifdef EOS_FORCE_2T
  integer,parameter :: overallCMask(EOSCOMP_NUM_COMPONENTS) = (/1,1,0/)
#else
  integer,parameter :: overallCMask(EOSCOMP_NUM_COMPONENTS) = (/1,1,1/)
#endif
  real,parameter :: eos_gammaRad = (4./3.)

  integer :: componentMask(EOSCOMP_NUM_COMPONENTS)

!!$  real,dimension(vecLen) :: gamcRow
  real,dimension(vecLen) :: gamIon, gamM1Ion, ggprodIon,ggprodInvIon,gam1InvIon
  real,dimension(vecLen) :: gamCEle, gamCM1Ele
  real,dimension(vecLen) :: gamEComb
  real :: ggprodEle
  real :: ggprodinvEle
  real :: gam1invEle, gam1invRad
  real :: dynamicZ, relA, Zinv, Zp, ZpInv
  real :: kBoltz
  real,dimension(NSPECIES) :: weight
  real :: rt,abarValue, abarInv, zbarValue, zbarFrac
  integer :: specieStart, specieEnd
  integer :: dens, temp, pres, eint, abar, zbar
  integer :: entr, dst, dsd
  integer :: dpt, dpd, det, ded, c_v, c_p, gamc, pel, ne, eta
  integer :: tempIon,tempEle,tempRad                                                
  integer :: eintIon,eintEle,eintRad                                                
  integer :: presIon,presEle,presRad
  integer :: entrEle, entrRad
  integer :: i, ilo,ihi
  integer :: seleVariant
  logical :: arrayBoundHackMode
  logical :: useNRForAddtlCombined !use vectors set by NR call for addtl. combined components
  !addtl. comp. = cp, cv, ...
  logical :: useNRForGamcCombined  !use vectors set by NR call for gamc combined components
  logical :: doUncoupledForAddtl ! assume components thermally uncoupled for addtl. out variables
  logical :: doUncoupledForGamc  ! assume components thermally uncoupled for gamc out variables
  logical :: doCoupledForGamc    ! handle components as thermally coupled for gamc out variables
  logical :: cvion_mask
  
  if(.not.present(massFrac)) then
     call Driver_abortFlash("[Eos] Multitype EOS implementation needs mass fractions")
  end if

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

  maskPlus(:) = maskPtr(:)
  maskPlus(EOS_PRESION) = .true.
  maskPlus(EOS_PRESELE) = .true.

  if (present(diagFlag)) diagFlag = 0

#ifdef DEBUG_EOS
  print*,'Calling Eos, mode',mode
#endif

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
!!$  rowLen = ihi - ilo + 1
#ifdef DEBUG_EOS
  if (ilo < 1 .OR. ilo > vecLen) then
     print*,'[Multitype/Eos] ilo is',ilo
     call Driver_abortFlash("[Multitype/Eos] invalid ilo")
  end if
  if (ihi < 1 .OR. ihi > vecLen) then
     print*,'[Multitype/Eos] ihi is',ihi
     call Driver_abortFlash("[Multitype/Eos] invalid ihi")
  end if
!!$  if (rowLen < 0 .OR. rowLen > vecLen) then
!!$     print*,'[Multitype/Eos] rowLen is',rowLen
!!$     call Driver_abortFlash("[Multitype/Eos] invalid rowLen")
!!$  end if
#endif
!!$  if (rowLen == 0) then
!!$     print*,'[Multitype/Eos] rowLen is 0.'
!!$  end if

  useNRForAddtlCombined = .FALSE.; useNRForGamcCombined = .FALSE.
  doUncoupledForAddtl = .FALSE.; doUncoupledForGamc = .FALSE.; doCoupledForGamc = .FALSE.

  kBoltz = eos_kerg
  seleVariant = 999 ! This effectively disables the setEleEntropy code below. - KW


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
  
#define dynamicZ eosData(zbar+ilo:zbar+ihi)
#define Zp (dynamicZ+1)
#define Zinv (1.0/dynamicZ)
#define ZpInv (1.0/Zp)


  select case (mode)
  case(MODE_DENS_TEMP_ALL)
  case(MODE_DENS_TEMP_COMP,MODE_DENS_TEMP_GATHER)
     useNRForAddtlCombined = .FALSE.; useNRForGamcCombined = .FALSE.
  case(MODE_DENS_TEMP,MODE_DENS_TEMP_EQUI)
!!$     useNRForAddtlCombined = .TRUE.; useNRForGamcCombined = .TRUE. !NONO in this file!
     doCoupledForGamc = .TRUE.

  case(MODE_DENS_EI_ALL)
  case(MODE_DENS_EI_COMP, MODE_DENS_EI_GATHER, &
       MODE_DENS_EI_RECAL_GATHER, &
       MODE_DENS_EI_SELE_GATHER)
     useNRForAddtlCombined = .FALSE.; useNRForGamcCombined = .FALSE.
  case(MODE_DENS_EI, MODE_DENS_EI_SCATTER)
     useNRForAddtlCombined = .TRUE.; useNRForGamcCombined = .TRUE.
     doCoupledForGamc = .TRUE.

    case(MODE_DENS_PRES_ALL)
    case(MODE_DENS_PRES_COMP)
     useNRForAddtlCombined = .FALSE.; useNRForGamcCombined = .FALSE.
    case(MODE_DENS_PRES)
     useNRForAddtlCombined = .TRUE.; useNRForGamcCombined = .TRUE.
     doCoupledForGamc = .TRUE.
  end select

  ! Some computations that are common to all cases handled below:
  ! Mostly, compute Abar and Zbar.
  do i = ilo,ihi
     specieStart = (i-1)*NSPECIES + 1
     specieEnd = i*NSPECIES

     call Multispecies_getSumInv(A, abarInv ,massFrac(specieStart:specieEnd))
     abarValue = 1.e0 / abarInv

     call Multispecies_getSumFrac(Z,zbarFrac,massFrac(specieStart:specieEnd))
     zbarValue = abarValue*zbarFrac

#ifdef ORIGINAL_GAMC_AVERAGE
     weight = massFrac(specieStart:specieEnd)*eos_gammam1j
     call Multispecies_getSumInv(A,rt, weight)
#else
     weight = massFrac(specieStart:specieEnd)*eos_gc
     call Multispecies_getSumInv(A, gmc, weight)
#endif

     eosData(abar+i) = abarValue
     eosData(zbar+i) = zbarValue
     if (mode==MODE_DENS_EI_RECAL_GATHER) then
#ifdef EOS_FORCE_2T
        call Hydro_recalibrateEints(eosData(eint+i), &
             eosData(eintIon+i),eosData(eintEle+i))
#else
        call Hydro_recalibrateEints(eosData(eint+i), &
             eosData(eintIon+i),eosData(eintEle+i),eosData(eintRad+i))
#endif
     end if
#ifdef ORIGINAL_GAMC_AVERAGE
     gamIon(i) = 1.0e0 + 1.0e0/(rt*eosData(abar+i))
     eosData(gamc+i) = 1.0e0 + 1.0e0/(rt*eosData(abar+i))
#else
     gamIon(i) = gmc*eosData(abar+i)
#endif
     gamM1Ion(i) = 1./(gamIon(i)-1.0)
  end do

  eosData(gamc+ilo:gamc+ihi) = (gamIon(ilo:ihi)+dynamicZ*eos_gammaEle)*ZpInv

#ifdef DEBUG_EOS
  print*,'eosdata(gamc...) at top of Eos.F90:',eosData(gamc+ilo:gamc+vecLen)
#endif

  select case (mode)
  case(MODE_DENS_TEMP,MODE_DENS_TEMP_COMP,MODE_DENS_TEMP_EQUI, &
       MODE_DENS_TEMP_ALL,MODE_DENS_TEMP_GATHER)
     ggprodIon = gamM1Ion * eos_gasConstant
     ggprodEle = eos_gammam1Ele * eos_gasConstant

  ! density, internal energies (or electron entropy) taken as input:
  case(MODE_DENS_EI, MODE_DENS_EI_SCATTER, &
       MODE_DENS_EI_ELE, &
       MODE_DENS_EI_COMP, MODE_DENS_EI_GATHER, &
       MODE_DENS_EI_RECAL_GATHER, &
       MODE_DENS_EI_ALL, &
       MODE_DENS_EI_SELE_GATHER)
     ggprodIon = gamM1Ion * eos_gasConstant
     ggprodEle = eos_gammam1Ele * eos_gasConstant

     ggprodinvIon = 1. / ggprodIon
     gam1invIon   = 1. / gamM1Ion
     ggprodinvEle = 1. / ggprodEle
     gam1invEle   = 1. / eos_gammam1Ele
!     gam1invRad   = 1. / eos_gammam1Rad


  ! density, pressures taken as input:
  case (MODE_DENS_PRES)
  case default
     call Driver_abortFlash("[Eos] Unrecognized input mode given to Eos")
  end select


     gamEComb(ilo:ihi) = (gamM1Ion(ilo:ihi)+dynamicZ*eos_gammam1Ele)*ZpInv
     gamEComb(ilo:ihi) = 1.0 + 1.0/gamEComb(ilo:ihi)

  if (doCoupledForGamc) then
     eosData(gamc+ilo:gamc+ihi) = gamEComb(ilo:ihi)
  end if

!============================================================================

  !! For allocatable arrays, set them up now.
#ifndef FIXEDBLOCKSIZE
  call eos_vecAlloc(vecLen)
#endif


!!$  print 9999,'Ye = ',Ye
!!$  9999 format(A7,1PG25.19)

#define Ye_ARRAYREFERENCE (eosData(zbar+ilo:zbar+ihi)/eosData(abar+ilo:abar+ihi))
#define ggprod (eos_gasConstant/(eosData(gamc+ilo:gamc+ihi)-1.0))

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


     call eos_multiTypeByTemp(MODE_DENS_TEMP,vecLen,eosData,ilo,ihi,massFrac,mask=maskPtr,componentMask=overallCMask)


  case (MODE_DENS_TEMP_COMP,MODE_DENS_TEMP_GATHER,MODE_DENS_TEMP_ALL)
     useNRForAddtlCombined = .FALSE.; useNRForGamcCombined = .FALSE.
     doUncoupledForAddtl = .TRUE.; doUncoupledForGamc = .TRUE.

     call eos_multiTypeByTemp(MODE_DENS_TEMP_ION,vecLen,eosData,ilo,ihi,massFrac,mask=maskPlus,componentMask=(/1,0,0/))
     gamCEle(ilo:ihi) = eosData(presIon+ilo:presIon+ihi)
     gamIon(ilo:ihi) = eosData(gamc+ilo:gamc+ihi)

     ! Work around. This is needed because EOS multiTypeByTemp will overwrite the
     ! the ion specific heat with 0.0 when called with MODE_DENS_TEMP_ELE
     cvion_mask = maskPlus(EOS_CVION)
     maskPlus(EOS_CVION) = .FALSE.

     call eos_multiTypeByTemp(MODE_DENS_TEMP_ELE,vecLen,eosData,ilo,ihi,massFrac,mask=maskPlus,componentMask=(/0,1,0/))
     eosData(presIon+ilo:presIon+ihi) = gamCEle(ilo:ihi)
     gamCEle(ilo:ihi) = eosData(gamc+ilo:gamc+ihi)

     maskPlus(EOS_CVION) = cvion_mask

     eosData(presRad+ilo:presRad+ihi)    = eos_ao3 * eosData(tempRad+ilo:tempRad+ihi)**4

     eosData(eintRad+ilo:eintRad+ihi)    = 3.0e0 * eosData(presRad+ilo:presRad+ihi) / eosData(dens+ilo:dens+ihi)

     eosData(entr+ilo:entr+ihi) = 0.0
     if (mode==MODE_DENS_TEMP_GATHER) then
!!$        if (ANY(maskPtr((/EOS_PRESELE,EOS_PRESION,EOS_PRESRAD/)) .EQV. .FALSE.)) then
!!$           call Driver_abortFlash("[Eos] cannot calculate MODE_DENS_TEMP_GATHER without component pressure masks.&
!!$                & Set mask appropriately.")
!!$        end if
        eosData(eint+ilo:eint+ihi) = ( &
             eosData(eintIon+ilo:eintIon+ihi)+ &
             eosData(eintEle+ilo:eintEle+ihi)+ &
             eosData(eintRad+ilo:eintRad+ihi))

        if(eos_combinedTempRule==0) then
           call eos_newtonRaphson(vecLen, eosData, ilo,ihi, MODE_DENS_EI, WANT_ENER, massFrac, maskPlus,cMask=overallCMask, &
                forceConstantInput=.TRUE.)
        endif

        ! Next lines moved from end of eos_newtonRaphson; pres changed - KW 2010-12-13
        eosData(pres+ilo:pres+ihi) = ( &
             eosData(presIon+ilo:presIon+ihi)+ &
             eosData(presEle+ilo:presEle+ihi)+ &
             eosData(presRad+ilo:presRad+ihi))
#ifdef EOS_FORCE_2T
        eosData(gamc+ilo:gamc+ihi)=eos_gammaEle
#else
        call setCombinedGamc(gamc,coupled=.FALSE.)
!        eosData(gamc+ilo:gamc+ihi)=gamcRow(ilo:ihi) !replace with expression based on component states!
#endif
        call setCombinedTemp()
!!!$!!!        eosData(entr+ilo:entr+ihi)=stotRow(ilo:ihi) !(replace with expression based on component states)

     else if (mode==MODE_DENS_TEMP_ALL) then
        doUncoupledForAddtl = .FALSE.; doUncoupledForGamc = .FALSE.

        call eos_multiTypeByTemp(MODE_DENS_TEMP,vecLen,eosData,ilo,ihi,massFrac,mask=maskPtr,componentMask=overallCMask)
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

#ifdef INCLUDE_INAPPROPRIATELY_IDEAL_CODE
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
#endif

     if(maskPtr(EOS_ENTRELE)) &
          call setEleEntropy(entrele,temp)

     call eos_multiTypeByTemp(MODE_DENS_TEMP_EQUI,vecLen,eosData,ilo,ihi,massFrac,mask=maskPtr,componentMask=overallCMask)

     eosData(entr+ilo:entr+ihi) = 0.0

     if(maskPtr(EOS_ENTRRAD)) &
          call setRadEntropy(entrrad,temp,eintRad)


#if(0)
  case (MODE_DENS_TEMP_ALL)
#ifdef INCLUDE_INAPPROPRIATELY_IDEAL_CODE
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
#endif

     if(maskPtr(EOS_ENTRELE)) &
          call setEleEntropy(entrele,tempEle)

     call eos_multiTypeByTemp(MODE_DENS_TEMP_ALL,vecLen,eosData,ilo,ihi,massFrac,mask=maskPtr,componentMask=overallCMask)

     eosData(entr+ilo:entr+ihi) = 0.0

     if(maskPtr(EOS_ENTRRAD)) &
          call setRadEntropy(entrrad,tempRad,eintRad)
#endif

  ! density, internal energy taken as input

  case (MODE_DENS_EI, MODE_DENS_EI_SCATTER)
     useNRForAddtlCombined = .TRUE.; useNRForGamcCombined = .TRUE.
     doUncoupledForAddtl = .FALSE.; doUncoupledForGamc = .FALSE.
     eosData(entr+ilo:entr+ihi) = 0.0                                 !overridden below!

     call eos_newtonRaphson(vecLen, eosData, ilo,ihi, mode, WANT_ENER, massFrac, maskPlus,cMask=overallCMask)
     
     ! Next lines moved from end of eos_newtonRaphson
!!$     eosData(pres+ilo:pres+ihi)=ptotRow(ilo:ihi)
!     eosData(gamc+ilo:gamc+ihi)=gamcRow(ilo:ihi)  ! Set in prologue! / set in epilogue!
     eosData(entr+ilo:entr+ihi)=stotRow(ilo:ihi)

     if (mode == MODE_DENS_EI_SCATTER) then
#if(0)
        ! should not be needed here, because done in eos_n.R. - was using tempRow
        eosData(tempIon+ilo:tempIon+ihi)=eosData(temp+ilo:temp+ihi)
        eosData(tempEle+ilo:tempEle+ihi)=eosData(temp+ilo:temp+ihi)
        eosData(tempRad+ilo:tempRad+ihi)=eosData(temp+ilo:temp+ihi)
#endif
!!$        if(maskPtr(EOS_PRESION)) eosData(presIon+ilo:presIon+ihi)=pCompRow(EOSCOMP_ION,ilo:ihi)
!!$        if(maskPtr(EOS_PRESELE)) eosData(presEle+ilo:presEle+ihi)=pCompRow(EOSCOMP_ELE,ilo:ihi)
!!$        if(maskPtr(EOS_PRESRAD)) eosData(presRad+ilo:presRad+ihi)=pCompRow(EOSCOMP_RAD,ilo:ihi)
!!$        if(maskPtr(EOS_EINTION)) eosData(eintIon+ilo:eintIon+ihi)=eCompRow(EOSCOMP_ION,ilo:ihi)
!!$        if(maskPtr(EOS_EINTELE)) eosData(eintEle+ilo:eintEle+ihi)=eCompRow(EOSCOMP_ELE,ilo:ihi)
!!$        if(maskPtr(EOS_EINTRAD)) eosData(eintRad+ilo:eintRad+ihi)=eCompRow(EOSCOMP_RAD,ilo:ihi)
        if(maskPtr(EOS_ENTRELE)) &
             call setEleEntropy(entrele,tempEle)
        if(maskPtr(EOS_ENTRRAD)) &
             call setRadEntropy(entrrad,tempRad,eintRad)
     else
        if(maskPtr(EOS_ENTRELE)) &
             call setEleEntropy(entrele,temp)
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
     call eos_newtonRaphson(vecLen, eosData, ilo,ihi, mode, WANT_ENER, massFrac, maskPlus,componentMask)



  case (MODE_DENS_EI_GATHER,MODE_DENS_EI_RECAL_GATHER)
     if (ANY(maskPtr((/EOS_PRESELE,EOS_PRESION,EOS_PRESRAD/)) .EQV. .FALSE.)) then
        call Driver_abortFlash("[Eos] cannot calculate MODE_DENS_EI_GATHER without component pressure masks.&
             & Set mask appropriately.")
     end if
     useNRForAddtlCombined = .FALSE.; useNRForGamcCombined = .FALSE.
     doUncoupledForAddtl = .TRUE.; doUncoupledForGamc = .TRUE. !or FALSE for _COMP?

     gamCEle(ilo:ihi) = eosData(eintRad+ilo:eintRad+ihi)
     maskPlus(EOS_ENTRELE) = .FALSE.
     maskPlus(EOS_ENTRRAD) = .FALSE.
     call eos_newtonRaphson(vecLen, eosData, ilo,ihi, MODE_DENS_EI_ION, WANT_ENER, massFrac, maskPlus,cMask=(/1,0,0/))
     gamIon(ilo:ihi) = eosData(gamc+ilo:gamc+ihi) !!  was gamcRow(ilo:ihi)
     maskPlus(EOS_ENTRELE) = maskPtr(EOS_ENTRELE)
     call eos_newtonRaphson(vecLen, eosData, ilo,ihi, MODE_DENS_EI_ELE, WANT_ENER, massFrac, maskPlus,cMask=(/0,1,0/))
     eosData(eintRad+ilo:eintRad+ihi) = gamCEle(ilo:ihi)
     gamCEle(ilo:ihi) = eosData(gamc+ilo:gamc+ihi) !! was gamcRow(ilo:ihi)

#ifdef DEBUG_MAY2011
!     print*,'gam1invIon e.o.Eos.F90:',gam1invIon
     print*,'presIon B Eos.F90:',eosData(presIon+1:presIon+vecLen)
     print*,'presEle B Eos.F90:',eosData(presEle+1:presEle+vecLen)
#endif

#ifdef EOS_FORCE_2T
     eosData(presRad+ilo:presRad+ihi) = 0.0
     eosData(tempRad+ilo:tempRad+ihi) = 0.0
#else
     eosData(presRad+ilo:presRad+ihi) = eosData(eintRad+ilo:eintRad+ihi) * eosData(dens+ilo:dens+ihi) / 3.0
     eosData(tempRad+ilo:tempRad+ihi) = sqrt(sqrt( eosData(presRad+ilo:presRad+ihi) / eos_ao3 ))
#endif

     eosData(eint+ilo:eint+ihi) = ( &
          eosData(eintIon+ilo:eintIon+ihi)+ &
          eosData(eintEle+ilo:eintEle+ihi)+ &
          eosData(eintRad+ilo:eintRad+ihi))

     eosData(entr+ilo:entr+ihi) = 0.0

     if(eos_combinedTempRule==0) then
        maskPlus(EOS_ENTRRAD) = maskPtr(EOS_ENTRRAD)
        call eos_newtonRaphson(vecLen, eosData, ilo,ihi, MODE_DENS_EI, WANT_ENER, massFrac, maskPlus,cMask=overallCMask)
        eosData(entr+ilo:entr+ihi)=stotRow(ilo:ihi) !(replace with expression based on component states)
     end if

     eosData(pres+ilo:pres+ihi) = ( &
          eosData(presIon+ilo:presIon+ihi)+ &
          eosData(presEle+ilo:presEle+ihi)+ &
          eosData(presRad+ilo:presRad+ihi))
     call setCombinedGamc(gamc,coupled=.FALSE.)

     call setCombinedTemp()

    !! Note that we make here a choice for the arbitrary additive constant in the electron entropy 
     if(maskPtr(EOS_ENTRELE)) &
          call setEleEntropy(entrele,tempEle)
     if(maskPtr(EOS_ENTRRAD)) &
          call setRadEntropy(entrrad,tempRad,eintRad)

#ifdef DEBUG_MAY2011
     print*,'presIon C Eos.F90:',eosData(presIon+ilo:presIon+ihi)
     print*,'presEle C Eos.F90:',eosData(presEle+ilo:presEle+ihi)
#endif



  case (MODE_DENS_EI_COMP)
     useNRForAddtlCombined = .FALSE.; useNRForGamcCombined = .FALSE.
     doUncoupledForAddtl = .TRUE.; doUncoupledForGamc = .FALSE.

     gamCEle(ilo:ihi) = eosData(eintRad+ilo:eintRad+ihi)
     call eos_newtonRaphson(vecLen, eosData, ilo,ihi, MODE_DENS_EI_ION, WANT_ENER, massFrac, maskPlus,cMask=(/1,0,0/))
     gamIon(ilo:ihi) = eosData(gamc+ilo:gamc+ihi) !!  was gamcRow(ilo:ihi)
     call eos_newtonRaphson(vecLen, eosData, ilo,ihi, MODE_DENS_EI_ELE, WANT_ENER, massFrac, maskPlus,cMask=(/0,1,0/))
     eosData(eintRad+ilo:eintRad+ihi) = gamCEle(ilo:ihi)
     gamCEle(ilo:ihi) = eosData(gamc+ilo:gamc+ihi) !! was gamcRow(ilo:ihi)

#ifdef EOS_FORCE_2T
     eosData(presRad+ilo:presRad+ihi) = 0.0
     eosData(tempRad+ilo:tempRad+ihi) = 0.0
#else
     eosData(presRad+ilo:presRad+ihi) = eosData(eintRad+ilo:eintRad+ihi) * eosData(dens+ilo:dens+ihi) / 3.0
     eosData(tempRad+ilo:tempRad+ihi) = sqrt(sqrt( eosData(presRad+ilo:presRad+ihi) / eos_ao3 ))
#endif

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

     call setCombinedGamc(gamc,coupled=.FALSE.,useCombPres=.FALSE.)

    !! Note that we make here a choice for the arbitrary additive constant in the electron entropy 
     if(maskPtr(EOS_ENTRELE)) &
          call setEleEntropy(entrele,tempEle)
     if(maskPtr(EOS_ENTRRAD)) &
          call setRadEntropy(entrrad,tempRad,eintRad)


  case (MODE_DENS_EI_ALL)    !DEV: inappropriately ideal for components.
     useNRForAddtlCombined = .TRUE.; useNRForGamcCombined = .TRUE. ! should DENS_EI_ALL work this way?
     doUncoupledForAddtl = .FALSE.; doUncoupledForGamc = .FALSE.
     eosData(presIon+ilo:presIon+ihi) = (1 * eosData(dens+ilo:dens+ihi)) * &
          eosData(eintIon+ilo:eintIon+ihi) * gam1invIon
     eosData(presEle+ilo:presEle+ihi) = (1 * eosData(dens+ilo:dens+ihi)) * &
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


     if (.TRUE.) then
        call eos_newtonRaphson(vecLen, eosData, ilo,ihi, MODE_DENS_EI, WANT_ENER, massFrac, maskPlus,cMask=overallCMask)
        ! Next lines moved from end of eos_newtonRaphson
!!$        eosData(pres+ilo:pres+ihi)=ptotRow(ilo:ihi) !Ok
!!$        eosData(gamc+ilo:gamc+ihi)=gamcRow(ilo:ihi) !Ok, I guess - if DENS_EI_ALL should work this way? - KW
        eosData(entr+ilo:entr+ihi)=stotRow(ilo:ihi) !(Ok, I guess - if DENS_EI_ALL should work this way? - KW)
     end if

    !! Note that we make here a choice for the arbitrary additive constant in the electron entropy 
     if(maskPtr(EOS_ENTRELE)) &
          call setEleEntropy(entrele,tempEle)
     if(maskPtr(EOS_ENTRRAD)) &
          call setRadEntropy(entrrad,tempRad,eintRad)


!!$  case (MODE_DENS_EI_ALL)    !inappropriately ideal.
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

     call eos_newtonRaphson(vecLen, eosData, ilo,ihi, MODE_DENS_PRES, WANT_PRES, massFrac, maskPlus, cMask=overallCMask)



  case (MODE_DENS_PRES_ION)     !all kinds of unwanted side effects... - KW
     call eos_newtonRaphson(vecLen, eosData, ilo,ihi, MODE_DENS_PRES, WANT_PRES, massFrac, maskPlus,cMask=(/1,0,0/))
!!$     eosData(tempIon+ilo:tempIon+ihi) = tempRow(ilo:ihi)
     eosData(tempIon+ilo:tempIon+ihi)=eosData(temp+ilo:temp+ihi)

  case (MODE_DENS_PRES_ELE)     !all kinds of unwanted side effects... - KW
     call eos_newtonRaphson(vecLen, eosData, ilo,ihi, MODE_DENS_PRES, WANT_PRES, massFrac, maskPlus,cMask=(/0,1,0/))
!!$     eosData(tempEle+ilo:tempEle+ihi) = tempRow(ilo:ihi)
     eosData(tempEle+ilo:tempEle+ihi)=eosData(temp+ilo:temp+ihi)

  case (MODE_DENS_PRES_RAD)
     eosData(tempRad+ilo:tempRad+ihi) = sqrt(sqrt( eosData(presRad+ilo:presRad+ihi) / eos_ao3 ))
     if(maskPtr(EOS_EINTRAD)) &
          eosData(eintRad+ilo:eintRad+ihi) = 3.0e0 * eosData(presRad+ilo:presRad+ihi) / eosData(dens+ilo:dens+ihi)


  case (MODE_DENS_PRES_COMP)     !all kinds of unwanted side effects... - KW
     call eos_newtonRaphson(vecLen, eosData, ilo,ihi, MODE_DENS_PRES, WANT_PRES, massFrac, maskPlus,cMask= (/1,0,0/))
!!$     eosData(tempIon+ilo:tempIon+ihi) = tempRow(ilo:ihi)
     eosData(tempIon+ilo:tempIon+ihi)=eosData(temp+ilo:temp+ihi)
     if(maskPtr(EOS_EINTION)) eosData(eintIon+ilo:eintIon+ihi)=eCompRow(EOSCOMP_ION,ilo:ihi)

     call eos_newtonRaphson(vecLen, eosData, ilo,ihi, MODE_DENS_PRES, WANT_PRES, massFrac, maskPlus,cMask= (/0,1,0/))
!!$     eosData(tempEle+ilo:tempEle+ihi) = tempRow(ilo:ihi)
     eosData(tempEle+ilo:tempEle+ihi)=eosData(temp+ilo:temp+ihi)
     if(maskPtr(EOS_EINTELE)) eosData(eintEle+ilo:eintEle+ihi)=eCompRow(EOSCOMP_ELE,ilo:ihi)

     eosData(tempRad+ilo:tempRad+ihi) = sqrt(sqrt( eosData(presRad+ilo:presRad+ihi) / eos_ao3 ))
     eosData(eintRad+ilo:eintRad+ihi)    = 3.0e0 * eosData(presRad+ilo:presRad+ihi) / eosData(dens+ilo:dens+ihi)

!!!!  case (MODE_DENS_PRES_ALL)

  case(MODE_DENS_EI_SELE_GATHER)
#ifdef EOS_FORCE_2T
     eosData(presRad+ilo:presRad+ihi) = 0.0
     eosData(tempRad+ilo:tempRad+ihi) = 0.0
#else
     ! We believe eintRad
     eosData(presRad+ilo:presRad+ihi) = eosData(eintRad+ilo:eintRad+ihi) * eosData(dens+ilo:dens+ihi) / 3.0
     eosData(tempRad+ilo:tempRad+ihi) = sqrt(sqrt( eosData(presRad+ilo:presRad+ihi) / eos_ao3 ))
#endif
     ! Done with _RAD (except entropy)

     ! We believe entrEle, do not believe eintEle
     gamIon(ilo:ihi) = eosData(eintRad+ilo:eintRad+ihi)
     call eos_newtonRaphson(vecLen, eosData, ilo,ihi, MODE_DENS_ENTR_ELE, WANT_ENTR, massFrac, maskPlus,cMask= (/0,1,0/))
#if(0)
     eosData(tempEle+ilo:tempEle+ihi) = exp( gam1invEle * ( eosData(entrele+ilo:entrele+ihi)/kBoltz &
                                                            + log( eosData(dens+ilo:dens+ihi) ) ) )
     eosData(eintEle+ilo:eintEle+ihi) = ggprodEle * eosData(tempEle+ilo:tempEle+ihi) &
          / eosData(abar+ilo:abar+ihi) * dynamicZ
     eosData(presEle+ilo:presEle+ihi) = eos_gasConstant * &
          eosData(dens+ilo:dens+ihi) * &
          eosData(tempEle+ilo:tempEle+ihi) / eosData(abar+ilo:abar+ihi) * dynamicZ
#endif
     gamCM1Ele(ilo:ihi) = eosData(presEle+ilo:presEle+ihi)
     gamCEle(ilo:ihi) = eosData(gamc+ilo:gamc+ihi) !! was gamcRow(ilo:ihi)
     ! eintEle and other _ELE have now been recomputed from entrEle

     if (ANY(maskPtr((/EOS_PRESELE,EOS_PRESION,EOS_PRESRAD/)) .EQV. .FALSE.)) then
        call Driver_abortFlash("[Eos] cannot calculate MODE_DENS_EI_GATHER without component pressure masks.&
             & Set mask appropriately.")
     end if

     ! We believe eint (combined), do NOT believe eintIon
     eosData(eintIon+ilo:eintIon+ihi) = eosData(eint+ilo:eint+ihi) - &
          eosData(eintEle+ilo:eintEle+ihi) - &
          eosData(eintRad+ilo:eintRad+ihi)
     ! eintIon has now been recomputed

     call eos_newtonRaphson(vecLen, eosData, ilo,ihi, MODE_DENS_EI_ION, WANT_ENER, massFrac, maskPlus,cMask=(/1,0,0/))
     eosData(eintRad+ilo:eintRad+ihi) = gamIon(ilo:ihi)
     gamIon(ilo:ihi) = eosData(gamc+ilo:gamc+ihi) !!  was gamcRow(ilo:ihi)
#if(0)
     eosData(presIon+ilo:presIon+ihi) = (1 * eosData(dens+ilo:dens+ihi)) * &
          eosData(eintIon+ilo:eintIon+ihi) * gam1invIon
     eosData(tempIon+ilo:tempIon+ihi) = eosData(eintIon+ilo:eintIon+ihi) * ggprodinvIon * &
          eosData(abar+ilo:abar+ihi)
#endif
     ! Other _ION have now been recomputed from eintIon

     eosData(presEle+ilo:presEle+ihi) = gamCM1Ele(ilo:ihi)

#if(0)
     eosData(eint+ilo:eint+ihi) = ( &
          eosData(eintIon+ilo:eintIon+ihi)+ &
          eosData(eintEle+ilo:eintEle+ihi)+ &
          eosData(eintRad+ilo:eintRad+ihi))
#endif
     useNRForAddtlCombined = .FALSE.; useNRForGamcCombined = .FALSE.
     doUncoupledForAddtl = .TRUE.; doUncoupledForGamc = .TRUE.

     ! Do Newton-Raphson iteration, which calls eos_multiTypeByTemp repeatedly.
     ! In this case, we call it mostly for the effect of computing an
     ! effective combined temperature for TEMP_VAR.  PRES_VAR also used to be
     ! set by this call.
     if (eos_combinedTempRule==0) then
        call eos_newtonRaphson(vecLen, eosData, ilo,ihi, MODE_DENS_EI, WANT_ENER, massFrac, maskPlus,cMask=overallCMask)
     end if
     
     ! Next lines moved from end of eos_newtonRaphson; pres changed - KW 2010-12-13
     eosData(pres+ilo:pres+ihi) = ( &
          eosData(presIon+ilo:presIon+ihi)+ &
          eosData(presEle+ilo:presEle+ihi)+ &
          eosData(presRad+ilo:presRad+ihi))
#ifdef EOS_FORCE_2T
     eosData(gamc+ilo:gamc+ihi)=eos_gammaEle
#else
     call setCombinedGamc(gamc,coupled=.TRUE.) ! Will be coupled=.FALSE., actually
!     eosData(gamc+ilo:gamc+ihi)=gamcRow(ilo:ihi) !replace with expression based on component states!
#endif
     call setCombinedTemp()
     eosData(entr+ilo:entr+ihi)=stotRow(ilo:ihi) !(replace with expression based on component states)

     if(maskPtr(EOS_ENTRRAD)) &
          call setRadEntropy(entrrad,tempRad,eintRad)

  ! unrecognized value for mode
  case default
     if (eos_meshMe==MASTER_PE) print*,"[Eos] Unrecognized input mode",mode
     call Logfile_stamp(mode,"[Eos] Unrecognized input mode given to Eos")
     call Driver_abortFlash("[Eos] Unrecognized input mode given to Eos")
  end select

!!$  print*,'eosData(abar+...) A.E.O. Eos():',eosData(abar+ilo:abar+ihi)
  if(useNRForGamcCombined) then
!!$     eosData(gamc+ilo:gamc+ihi) = gamcRow(ilo:ihi)
  else if (doUncoupledForGamc) then
     call setCombinedGamc(gamc,coupled=.FALSE.)
#if(0)
  else if (doUncoupledForGamc) then
     call averageV1Sc2ByPressureFraction(gamc, &
          gamIon, eos_gammaEle, eos_gammaRad)
#endif
  else
     ! leave as set above - KW
  end if

  if(present(mask)) then

     if(mask(EOS_DPT)) then
        dpt = (EOS_DPT-1)*vecLen
        if(useNRForAddtlCombined) then
           ! done in NR... !eosData(dpt+ilo:dpt+ihi) = dptRow(ilo:ihi)
        else if (doUncoupledForAddtl) then
           eosData(dpt+ilo:dpt+ihi) = &
                (  eos_gasConstant * eosData(dens+ilo:dens+ihi) * &
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
!!$     if(mask(EOS_DET))then
!!$        det = (EOS_DET-1)*vecLen
!!$        if(useNRForAddtlCombined) then
!!$           ! done in NR... !eosData(det+ilo:det+ihi) = detRow(ilo:ihi)
!!$           ! done in eos_multiTypeByTemp calls otherwise
!!$        end if
!!$     end if
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
           call Driver_abortFlash("[Eos] Cannot calculate EOS_DST without EOS_DET and EOS_DPT")
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
           call Driver_abortFlash("[Eos] Cannot calculate EOS_DSD without EOS_DED and EOS_DPD")
        end if
     end if


     if(mask(EOS_PEL))then 
        pel = (EOS_PEL-1)*vecLen
!!$        eosData(pel+ilo:pel+ihi) = 0.
     end if
     if(mask(EOS_NE))then 
        ne = (EOS_NE-1)*vecLen
!!$        eosData(ne+ilo:ne+ihi) = 0.
     end if
     if(mask(EOS_ETA))then 
        call Driver_abortFlash("[Eos] cannot calculate ETA in the multiTemp / Gamma implementation.")
        eta = (EOS_ETA-1)*vecLen
        eosData(eta+ilo:eta+ihi) = 0.
     end if
     
     if(mask(EOS_CV))then
        if(mask(EOS_DET)) then
           c_v = (EOS_CV-1)*vecLen
           det = (EOS_DET-1)*vecLen
           if(useNRForAddtlCombined) then
              eosData(c_v+ilo:c_v+ihi) = cvRow(ilo:ihi)
           else if (doUncoupledForAddtl) then
              eosData(c_v+ilo:c_v+ihi) = eosData(det+ilo:det+ihi)
           else
              eosData(c_v+ilo:c_v+ihi) = eosData(det+ilo:det+ihi)
           end if
        else
           call Driver_abortFlash("[Eos] cannot calculate C_V without DET.  Set mask appropriately.")
        end if
     end if
     ! ideal gas -- all gammas are equal
     if(mask(EOS_CP))then
        if(mask(EOS_CV).and.mask(EOS_DET)) then
           c_p = (EOS_CP-1)*vecLen
           if(useNRForAddtlCombined) then
              eosData(c_p+ilo:c_p+ihi) = cpRow(ilo:ihi)
           else if (doUncoupledForAddtl) then
              call averageV1Sc2ByPressureFraction(c_p, &
                   gamIon, eos_gammaEle, eos_gammaRad)
              eosData(c_p+ilo:c_p+ihi) = eosData(c_p+ilo:c_p+ihi)*eosData(c_v+ilo:c_v+ihi)
           else
              eosData(c_p+ilo:c_p+ihi) = eosData(gamc+ilo:gamc+ihi)*eosData(c_v+ilo:c_v+ihi)
           end if
        else
           call Driver_abortFlash("[Eos] cannot calculate C_P without C_V and DET.  Set mask appropriately.")
        end if
        
     end if

!      if(mask(EOS_CVION))then
!            c_v = (EOS_CVION-1)*vecLen
!            eosData(c_v+ilo:c_v+ihi) = ggprod  / eosData(abar+ilo:abar+ihi)
!      end if
!      if(mask(EOS_CVELE))then
! !!$           c_v = (EOS_CVELE-1)*vecLen
! !!$           eosData(c_v+ilo:c_v+ihi) = &
! !!$                ggprod * dynamicZ / eosData(abar+ilo:abar+ihi)
!      end if
  end if


  if (arrayBoundHackMode .eqv. .true.) then
     deallocate(maskPtr)
  end if


  !! Close up arrays if previously allocated
#ifndef FIXEDBLOCKSIZE  
  call eos_vecDealloc()
#endif

contains
  subroutine setEleEntropy(entrComponent,tempComponent)
  use eos_mgammaData, ONLY : eos_eMass
  use eos_helmConstData, ONLY: eos_avo,eos_kergavo, eos_h,eos_hbar
    integer,intent(IN) :: entrComponent, tempComponent

    select case (seleVariant)
    case(1)
    eosData(entrComponent+ilo:entrComponent+ihi) = &
         kBoltz * &
         (eos_gammam1Ele * log( eosData(tempComponent+ilo:tempComponent+ihi) ) &
                      -    log( eosData(dens+ilo:dens+ihi)) )
    case(2)
    eosData(entrComponent+ilo:entrComponent+ihi) = &
         eos_kergavo * eosData(zbar+ilo:zbar+ihi)/eosData(abar+ilo:abar+ihi) * kBoltz * &
         (eos_gammam1Ele * log( eosData(tempComponent+ilo:tempComponent+ihi) * &
                                kBoltz *eos_eMass / (eos_h*eos_hbar) ) &
                      -    log( eosData(dens+ilo:dens+ihi) * eos_avo *eosData(zbar+ilo:zbar+ihi) / eosData(abar+ilo:abar+ihi)) )
    case(3)
    eosData(entrComponent+ilo:entrComponent+ihi) = &
         eos_kergavo * Ye_ARRAYREFERENCE * kBoltz * &
         (eos_gammam1Ele * log( eosData(tempComponent+ilo:tempComponent+ihi) * &
                                kBoltz *eos_eMass / (eos_h*eos_hbar) ) &
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
         (                 log( eosData(tempComponent+ilo:tempComponent+ihi)**eos_gammam1Ele &
                           /    eosData(dens+ilo:dens+ihi)                    ) )
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
    use eos_idealGammaData, ONLY : eos_eMass
    use eos_helmConstData, ONLY: eos_avo,eos_kergavo, eos_h,eos_hbar
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

  subroutine setCombinedGamc(gamcCombined,coupled,useCombPres)
    integer, intent(in) :: gamcCombined
    logical, intent(in) :: coupled
    logical, intent(in),OPTIONAL :: useCombPres
    logical :: useCombP
    useCombP = .TRUE.
    if (present(useCombPres)) useCombP = useCombPres

    if (coupled) then
! DEV: eos_gammam1Ele -> gamIon(ilo:ihi), considerations?
       eosData(gamcCombined+ilo:gamcCombined+ihi) = &
         ( gamM1Ion(ilo:ihi)*eosData(presIon+ilo:presIon+ihi)+ &
           eos_gammam1Ele*eosData(presEle+ilo:presEle+ihi)+ &
           eos_gammam1Rad*eosData(presRad+ilo:presRad+ihi) ) / eosData(pres+ilo:pres+ihi)
       eosData(gamcCombined+ilo:gamcCombined+ihi) = &
            1.0 + 1.0/eosData(gamcCombined+ilo:gamcCombined+ihi) 
#ifdef DEBUG_EOS
       print*,'setCGamcC>',eosData(gamcCombined+1:gamcCombined+vecLen)
#endif
    else if (useCombP) then
       eosData(gamcCombined+ilo:gamcCombined+ihi) = &
         ( gamIon(ilo:ihi)*eosData(presIon+ilo:presIon+ihi)+ &
           gamCEle(ilo:ihi)*eosData(presEle+ilo:presEle+ihi)+ &
           eos_gammaRad*eosData(presRad+ilo:presRad+ihi) ) / eosData(pres+ilo:pres+ihi)
#ifdef DEBUG_EOS
       print*,'setCGamcU>',eosData(gamcCombined+1:gamcCombined+vecLen)
#endif
    else
       eosData(gamcCombined+ilo:gamcCombined+ihi) = &
         ( gamIon(ilo:ihi)*eosData(presIon+ilo:presIon+ihi)+ &
           gamCEle(ilo:ihi)*eosData(presEle+ilo:presEle+ihi)+ &
           eos_gammaRad*eosData(presRad+ilo:presRad+ihi) ) / &
	   (eosData(presIon+ilo:presIon+ihi)+eosData(presEle+ilo:presEle+ihi)+eosData(presRad+ilo:presRad+ihi))
#ifdef DEBUG_EOS
       print*,'setCGamcu>',eosData(gamcCombined+1:gamcCombined+vecLen)
#endif
    end if

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
    real, intent(in),dimension(1:vecLen) :: compIon,compEle,compRad

    eosData(combined+ilo:combined+ihi) = &
         ( compIon(ilo:ihi)*eosData(presIon+ilo:presIon+ihi)+ &
           compEle(ilo:ihi)*eosData(presEle+ilo:presEle+ihi)+ &
           compRad(ilo:ihi)*eosData(presRad+ilo:presRad+ihi) ) / eosData(pres+ilo:pres+ihi)

  end subroutine averageVectRByPressureFraction

  subroutine averageV1Sc2ByPressureFraction(combined,compIon,scal2,scal3)
    integer, intent(in) :: combined
    real, intent(in),dimension(1:vecLen) :: compIon
    real, intent(in) :: scal2,scal3

    eosData(combined+ilo:combined+ihi) = &
         ( compIon(ilo:ihi)*eosData(presIon+ilo:presIon+ihi)+ &
           scal2            *eosData(presEle+ilo:presEle+ihi)+ &
           scal3            *eosData(presRad+ilo:presRad+ihi) ) / eosData(pres+ilo:pres+ihi)

  end subroutine averageV1Sc2ByPressureFraction



end subroutine Eos

  subroutine eos_newtonRaphson(vecLen, eosData, ilo,ihi, mode, whatWanted, &
       massFrac, eosMaskPlus,cMask,forceConstantInput)
  use eos_mtInterface, ONLY : eos_multiTypeByTemp
  use Eos_data, ONLY : eos_smallT, eos_largeT, eos_tol, eos_maxNewton, &
       eos_forceConstantInput
  use Eos_mgammaData, ONLY : eos_maxFactorUp, eos_maxFactorDown
  use eos_vecData, ONLY:  denRow
    implicit none
!#define DEBUG_EOS
    integer,intent(IN) :: vecLen, mode, ilo,ihi
    real,INTENT(inout), dimension(EOS_NUM*vecLen) :: eosData 
    integer,intent(IN) :: whatWanted

    real,   INTENT(in),dimension(NSPECIES*vecLen) :: massFrac
    logical,optional, dimension(EOS_VARS+1:EOS_NUM),INTENT(INOUT)::eosMaskPlus
    integer,optional, dimension(EOSCOMP_NUM_COMPONENTS),INTENT(in)::cMask
    logical,optional,intent(in) :: forceConstantInput

    integer :: ancMode
    integer :: k
    integer :: vecBegin,vecEnd  !DEV: needed?
    logical :: forceConstInp
  ! declare some local storage for the temperature during Newton iterations
    real,dimension(vecLen)::  tempRow
  ! declare some local storage for the results of the Newton iteration
    real,dimension(vecLen)::  ewantRow, tnew, error,pwantRow, swantRow
    real,dimension(vecLen)::  etotRow, &
       tempRadRow,tempIonRow,tempEleRow
    real,dimension(vecLen)::  ptotRow
    real,dimension(vecLen)::  stotRow
    target :: ewantRow,pwantRow,swantRow
    target :: etotRow,ptotRow,stotRow
    integer :: dens, temp, pres, eint
    integer :: entr, dst
    integer :: dpt, det
    integer :: tempIon,tempEle,tempRad                                                
    integer :: eintIon,eintEle,eintRad
    integer :: presIon,presEle,presRad
    integer :: entrEle, entrRad
    integer :: i
  !  local storage for forcedConstantInput -- could be allocatable, but might be so slow
  !  that it's not worth the small storage save.
    real,dimension(vecLen)::  psaveRow, esaveRow

    integer :: spec_num
    character(len=MAX_STRING_LENGTH) :: spec_str
    real, pointer, dimension(:) :: xXtotRow, xXwantRow
    real    :: told
    real,dimension(vecLen):: bnd_lo, bnd_hi


    if (present(forceConstantInput)) then
       forceConstInp = forceConstantInput
    else
       forceConstInp = eos_forceConstantInput
    end if
    vecBegin = ilo
    vecEnd   = ihi

  pres = (EOS_PRES-1)*vecLen
  dens = (EOS_DENS-1)*vecLen
  temp = (EOS_TEMP-1)*vecLen
  eint = (EOS_EINT-1)*vecLen
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

    do k = ilo,ihi
       tempRow(k)    = eosData(temp+k)
       denRow(k)     = eosData(dens+k)
       ewantRow(k)   = eosData(eint+k)   ! store desired internal energy for MODE_DENS_EI etc. cases
       pwantRow(k)   = eosData(pres+k)   ! store desired pressure for MODE_DENS_PRES etc. cases
    end do



    select case (mode)
    case(MODE_DENS_TEMP_GATHER)
       ancMode = MODE_DENS_TEMP
    case(MODE_DENS_TEMP,MODE_DENS_TEMP_COMP,MODE_DENS_TEMP_EQUI, &
         MODE_DENS_TEMP_ALL)
       ancMode = -1
       ! density, internal energies (or electron entropy) taken as input:
    case(MODE_DENS_EI)
       ancMode = MODE_DENS_TEMP
       eosMaskPlus(EOS_PRESION) = .FALSE.
       eosMaskPlus(EOS_PRESELE) = .FALSE.
    case(MODE_DENS_EI_GATHER, &
         MODE_DENS_EI_RECAL_GATHER, &
         MODE_DENS_EI_ALL, &
         MODE_DENS_EI_SELE_GATHER)
       ancMode = MODE_DENS_TEMP

!!$    case(MODE_DENS_EI_ALL)
!!$       ancMode = MODE_DENS_TEMP_ALL

    case(MODE_DENS_EI_SCATTER)
       ancMode = MODE_DENS_TEMP_EQUI

    case(MODE_DENS_EI_ION)
       ancMode = MODE_DENS_TEMP_ION
    case(MODE_DENS_EI_ELE)
       ancMode = MODE_DENS_TEMP_ELE
    case(MODE_DENS_EI_RAD)
       ancMode = MODE_DENS_TEMP_RAD

    case(MODE_DENS_EI_COMP)
!!$       ancMode = MODE_DENS_TEMP_GATHER
       ancMode = -1

       ! density, pressures taken as input:
    case(MODE_DENS_PRES)
       ancMode = MODE_DENS_TEMP
       eosMaskPlus(EOS_PRESION) = .FALSE.
       eosMaskPlus(EOS_PRESELE) = .FALSE.
!!!       ancMode = MODE_DENS_TEMP_EQUI  ! Works, too, for PRES_TEMP!
!!!       eosMaskPlus(EOS_PRESION) = .FALSE.
!!!       eosMaskPlus(EOS_PRESELE) = .FALSE.

    case(MODE_DENS_PRES_ALL)
       ancMode = MODE_DENS_TEMP
!!!       ancMode = MODE_DENS_TEMP_EQUI  ! Works, too, for PRES_TEMP!

!!$    case(MODE_DENS_PRES_ALL)
!!$       ancMode = MODE_DENS_TEMP_ALL
!!$
!!$    case(MODE_DENS_PRES_SCATTER)
!!$       ancMode = MODE_DENS_TEMP_EQUI

    case(MODE_DENS_PRES_COMP)
       ancMode = MODE_DENS_TEMP_GATHER

    case(MODE_DENS_PRES_ION)
       ancMode = MODE_DENS_TEMP_ION
    case(MODE_DENS_PRES_ELE)
       ancMode = MODE_DENS_TEMP_ELE
    case(MODE_DENS_PRES_RAD)
       ancMode = MODE_DENS_TEMP_RAD

    case(MODE_DENS_ENTR_ELE)
       ancMode = MODE_DENS_TEMP_ELE

    case default
       ancMode = MODE_DENS_TEMP_ALL
    end select


    select case (mode)
    case(MODE_DENS_EI_ION)
       eosMaskPlus(EOS_EINTION) = .true.
!!$       eosMaskPlus(EOS_PRESELE) = .FALSE.  ! DEV: not sure why commented out in r14167 - KW
       eosMaskPlus(EOS_PRESRAD) = .FALSE.
       ewantRow(vecBegin:vecEnd) = eosData(eintIon+vecBegin:eintIon+vecEnd)
    case(MODE_DENS_EI_ELE)
       eosMaskPlus(EOS_EINTELE) = .true.
       eosMaskPlus(EOS_PRESION) = .FALSE.
       eosMaskPlus(EOS_PRESRAD) = .FALSE.
       ewantRow(vecBegin:vecEnd) = eosData(eintEle+vecBegin:eintEle+vecEnd)
    case(MODE_DENS_EI_RAD)
       eosMaskPlus(EOS_EINTRAD) = .true.
       eosMaskPlus(EOS_PRESION) = .FALSE.
       eosMaskPlus(EOS_PRESELE) = .FALSE.
       ewantRow(vecBegin:vecEnd) = eosData(eintRad+vecBegin:eintRad+vecEnd)
    case default   ! store desired internal energy for MODE_DENS_EI etc. cases
       ewantRow(vecBegin:vecEnd) = eosData(eint+vecBegin:eint+vecEnd)
    end select
    select case (mode)
    case(MODE_DENS_PRES_ION)
       eosMaskPlus(EOS_PRESION) = .true.
       pwantRow(vecBegin:vecEnd) = eosData(presIon+vecBegin:presIon+vecEnd)
    case(MODE_DENS_PRES_ELE)
       eosMaskPlus(EOS_PRESELE) = .true.
       pwantRow(vecBegin:vecEnd) = eosData(presEle+vecBegin:presEle+vecEnd)
    case(MODE_DENS_PRES_RAD)
       eosMaskPlus(EOS_PRESRAD) = .true.
       pwantRow(vecBegin:vecEnd) = eosData(presRad+vecBegin:presRad+vecEnd)
    case default
       pwantRow(vecBegin:vecEnd) = eosData(pres+vecBegin:pres+vecEnd)
    end select
    select case (mode)
    case(MODE_DENS_ENTR_ELE)
       eosMaskPlus(EOS_ENTRELE) = .true.
       swantRow(vecBegin:vecEnd) = eosData(entrEle+vecBegin:entrEle+vecEnd)
!       eosMaskPlus(EOS_PRESION) = .FALSE.
       eosMaskPlus(EOS_PRESRAD) = .FALSE.
    case(MODE_DENS_ENTR_RAD)
       eosMaskPlus(EOS_ENTRRAD) = .true.
       swantRow(vecBegin:vecEnd) = eosData(entrRad+vecBegin:entrRad+vecEnd)
    case default
       swantRow(vecBegin:vecEnd) = eosData(entr+vecBegin:entr+vecEnd)
    end select


    !! Save the input parameters if forceConstInp is requested
    if (forceConstInp) then
       esaveRow = ewantRow
       psaveRow = pwantRow
    end if

    if (whatWanted == WANT_ENER) then
       xXtotRow  => etotRow
       xXwantRow => ewantRow
       eosMaskPlus(EOS_DET) = .true.
    else if (whatWanted == WANT_ENTR) then
       xXtotRow  => ptotRow
       xXwantRow => pwantRow
       eosMaskPlus(EOS_DST) = .true.
       eosMaskPlus(EOS_DET) = .true.
       eosMaskPlus(EOS_DPT) = .true.
    else
       xXtotRow  => stotRow
       xXwantRow => swantRow
       eosMaskPlus(EOS_DPT) = .true.
    end if


    ! Initialize the errors
    error(:) = 0.0e0

    !! set wide bracket bounds
    bnd_lo(:) = eos_smallt
    bnd_hi(:) = eos_largeT

    ! Do the first eos call with all the zones in the pipe
    !  NOTE that eos_multiTypeByTemp can ONLY operate in the equivalent of
    !  MODE_DENS_TEMP, as it returns pressure, energy and derivatives only
    !  So if you send in a crappy temperature here, you'll get a crappy starting
    !  position and the iteration won't converge.
    !  Initial temperature here is what is stored in the grid, even though we 
    !    SUSPECT this is not in equilibrium (or we wouldn't be calling Eos if it was fine)
    select case (ancMode)
    case(MODE_DENS_TEMP_ION)
       tempRow(ilo:ihi) = eosdata(tempIon+ilo:tempIon+ihi)
    case(MODE_DENS_TEMP_ELE)
       tempRow(ilo:ihi) = eosdata(tempEle+ilo:tempEle+ihi)
    case(MODE_DENS_TEMP_RAD)
       tempRow(ilo:ihi) = eosdata(tempRad+ilo:tempRad+ihi)
    end select
    if (.NOT. all(tempRow(ilo:ihi) > 0)) &
         print*,'eos_nr WARN:',tempRow(ilo:ihi)

    do k = ilo, ihi
       if (tempRow(k) > eos_largeT) tempRow(k) = eos_largeT
       if (tempRow(k) < eos_smallT) tempRow(k) = eos_smallT
    enddo
    tempRow(ilo:ihi) = max(0.0,eos_smallT,tempRow(ilo:ihi))
    call local_setEosTempsFromTempRow(ilo,ihi,ancMode,cMask)
!!$    if (.NOT. present(cMask)) then
!!$       tempRadRow(ilo:ihi) = tempRow(ilo:ihi)
!!$       tempIonRow(ilo:ihi) = tempRow(ilo:ihi)
!!$       tempEleRow(ilo:ihi) = tempRow(ilo:ihi)
!!$    else
!!$       tempRadRow(ilo:ihi) = eosdata(tempRad+ilo:tempRad+ihi)
!!$       tempIonRow(ilo:ihi) = eosdata(tempIon+ilo:tempIon+ihi)
!!$       tempEleRow(ilo:ihi) = eosdata(tempEle+ilo:tempEle+ihi)
!!$       if (cMask(EOSCOMP_RAD).NE.0) tempRadRow(ilo:ihi) = tempRow(ilo:ihi)
!!$       if (cMask(EOSCOMP_ION).NE.0) tempIonRow(ilo:ihi) = tempRow(ilo:ihi)
!!$       if (cMask(EOSCOMP_ELE).NE.0) tempEleRow(ilo:ihi) = tempRow(ilo:ihi)
!!$    end if
!!$    select case (ancMode)
!!$    case(MODE_DENS_TEMP,MODE_DENS_TEMP_EQUI)
!!$       eosdata(temp+ilo:temp+ihi)       = tempRow(ilo:ihi)
!!$    case(MODE_DENS_TEMP_ALL)
!!$       eosdata(temp+ilo:temp+ihi)       = tempRow(ilo:ihi)
!!$       eosdata(tempIon+ilo:tempIon+ihi) = tempIonRow(ilo:ihi)
!!$       eosdata(tempEle+ilo:tempEle+ihi) = tempEleRow(ilo:ihi)
!!$       eosdata(tempRad+ilo:tempRad+ihi) = tempRadRow(ilo:ihi)
!!$    case(MODE_DENS_TEMP_GATHER)
!!$       eosdata(tempIon+ilo:tempIon+ihi) = tempIonRow(ilo:ihi)
!!$       eosdata(tempEle+ilo:tempEle+ihi) = tempEleRow(ilo:ihi)
!!$       eosdata(tempRad+ilo:tempRad+ihi) = tempRadRow(ilo:ihi)
!!$    case(MODE_DENS_TEMP_ION)
!!$       eosdata(tempIon+ilo:tempIon+ihi) = tempIonRow(ilo:ihi)
!!$    case(MODE_DENS_TEMP_ELE)
!!$       eosdata(tempEle+ilo:tempEle+ihi) = tempEleRow(ilo:ihi)
!!$    case(MODE_DENS_TEMP_RAD)
!!$       eosdata(tempRad+ilo:tempRad+ihi) = tempRadRow(ilo:ihi)
!!$    end select

!!$    call eos_multiTypeByTemp(vecBegin,vecEnd,mask=eosMaskPlus,componentMask=cMask)
    call eos_multiTypeByTemp(ancMode,vecLen,eosData,vecBegin,vecEnd,massFrac,mask=eosMaskPlus,componentMask=cMask)
    call local_setTempRow(vecBegin,vecEnd,ancMode)
!!$    select case (ancMode)
!!$       case(MODE_DENS_TEMP_ION)
!!$          tempRow(vecBegin:vecEnd) = eosData(tempIon+vecBegin:tempIon+vecEnd)
!!$       case(MODE_DENS_TEMP_ELE)
!!$          tempRow(vecBegin:vecEnd) = eosData(tempEle+vecBegin:tempEle+vecEnd)
!!$       case(MODE_DENS_TEMP_RAD)
!!$          tempRow(vecBegin:vecEnd) = eosData(tempRad+vecBegin:tempRad+vecEnd)
!!$       case default
!!$       tempRow(vecBegin:vecEnd) = eosData(temp+vecBegin:temp+vecEnd) !maybe not..
!!$    end select

    call local_setTots(vecBegin,vecEnd,mode)
!!$    select case (mode)
!!$       case(MODE_DENS_EI_ION,MODE_DENS_PRES_ION)
!!$          etotRow(vecBegin:vecEnd) = eosData(eintIon+vecBegin:eintIon+vecEnd)
!!$          ptotRow(vecBegin:vecEnd) = eosData(presIon+vecBegin:presIon+vecEnd)
!!$       case(MODE_DENS_EI_ELE,MODE_DENS_PRES_ELE,MODE_DENS_ENTR_ELE)
!!$          etotRow(vecBegin:vecEnd) = eosData(eintEle+vecBegin:eintEle+vecEnd)
!!$          ptotRow(vecBegin:vecEnd) = eosData(presEle+vecBegin:presEle+vecEnd)
!!$          stotRow(vecBegin:vecEnd) = eosData(entrEle+vecBegin:entrEle+vecEnd)
!!$       case(MODE_DENS_EI_RAD,MODE_DENS_PRES_RAD)
!!$          etotRow(vecBegin:vecEnd) = eosData(eintRad+vecBegin:eintRad+vecEnd)
!!$          ptotRow(vecBegin:vecEnd) = eosData(presRad+vecBegin:presRad+vecEnd)
!!$       case default
!!$          etotRow(vecBegin:vecEnd) = eosData(eint+vecBegin:eint+vecEnd)
!!$          ptotRow(vecBegin:vecEnd) = eosData(pres+vecBegin:pres+vecEnd)
!!$       end select

!!$    gamcRow(vecBegin:vecEnd) = eosData(gamc+vecBegin:gamc+vecEnd)
!!!$!!!    stotRow(vecBegin:vecEnd) = eosData(entr+vecBegin:entr+vecEnd)
!!$    print*,'ptotRow(NR1(Eos)):',ptotRow(vecBegin:vecEnd)
    !  Now eos_multiTypeByTemp has returned ptotRow, etotRow, dXXtRow, and gamcRow


    !  Create initial condition
    do k = vecBegin, vecEnd

       if (whatWanted == WANT_ENER) then
          det = (EOS_DET-1)*vecLen
          !  ewantRow is our desired EI input
          tnew(k) = tempRow(k) - (etotRow(k) - ewantRow(k)) / eosData(det+k)
       else if (whatWanted == WANT_ENTR) then
          dst = (EOS_DST-1)*vecLen
          !  swantRow is our desired entropy input
          tnew(k) = tempRow(k) - (stotRow(k) - swantRow(k)) / eosData(dst+k)
       else                    !pressure is wanted
          dpt = (EOS_DPT-1)*vecLen
          tnew(k) = tempRow(k) - (ptotRow(k) - pwantRow(k)) / eosData(dpt+k)
       end if

       ! if we jump out of the brackets, check if there is a solution
       if ( tnew(k) < bnd_lo(k) ) then
          told = tempRow(k)
          tempRow(k) = bnd_lo(k)
          call local_setEosTempsFromTempRow(k,k,ancMode,cMask)
          call eos_multiTypeByTemp(ancMode,vecLen,eosData,k,k,massFrac,mask=eosMaskPlus,componentMask=cMask)
          call local_setTempRow(k,k,ancMode)
          call local_setTots(k,k,mode)
          ! assume eint (or pres) is a monotonically increasing function of T
          if ( xXwantRow(k) <= xXtotRow(k) ) then
             ! this will force exit from the loop via error=0
             tnew(k) = tempRow(k)
          else
             ! bisection
             tempRow(k) = told
             bnd_hi(k) = told
             tnew(k) = 0.5*(bnd_lo(k)+bnd_hi(k))
          endif
       else if ( tnew(k) > bnd_hi(k) ) then
          told = tempRow(k)
          tempRow(k) = bnd_hi(k)
          call local_setEosTempsFromTempRow(k,k,ancMode,cMask)
          call eos_multiTypeByTemp(ancMode,vecLen,eosData,k,k,massFrac,mask=eosMaskPlus,componentMask=cMask)
          call local_setTempRow(k,k,ancMode)
          call local_setTots(k,k,mode)
          ! assume eint (or pres) is a monotonically increasing function of T
          if ( xXwantRow(k) >= xXtotRow(k) ) then
             ! this will force exit from the loop via error=0
             tnew(k) = tempRow(k)
          else
             ! bisection
             tempRow(k) = told
             bnd_lo(k) = told
             tnew(k) = 0.5*(bnd_lo(k)+bnd_hi(k))
          endif
       else
          ! using derivative as an indicator of whether the old guess was
          ! an upper lower bound, update bracket range
          ! (assumes derivative is always well-behaved)
          if ( tnew(k) >= tempRow(k) ) then
             bnd_lo(k) = tempRow(k)
          else
             bnd_hi(k) = tempRow(k)
          endif
       endif

       ! Don't allow the temperature to change by more than an order of magnitude 
       ! in a single iteration
#ifdef DEBUG_EOS
       print*,'MAX_Fs:',eos_maxFactorUp,eos_maxFactorDown
#endif
       if (tnew(k) .GT. eos_maxFactorUp*tempRow(k)) tnew(k) =  & 
            &           eos_maxFactorUp*tempRow(k)
       if (tnew(k) .LT. eos_maxFactorDown*tempRow(k)) tnew(k) =  & 
            &           eos_maxFactorDown*tempRow(k)

       ! Compute the error
       error(k) = abs((tnew(k) - tempRow(k)) / tempRow(k))

#ifdef DEBUG_EOS
       if (whatWanted == WANT_ENTR) then
          print*,'E0:',error(k),(stotRow(k) - swantRow(k)),k,tnew(k),tempRow(k) 
       else
       print*,'error  (NR0(Eos)):',error(k),k,tnew(k),tempRow(k)
       endif
#endif

       ! Store the new temperature
       tempRow(k) = tnew(k)

       ! DEV: Disable the following 'freezing' logic? - KW
       ! Check if we are freezing, if so set the temperature to smallt, and adjust 
       ! the error so we don't wait for this one
       if (tempRow(k) .LT. eos_smallt) then
          tempRow(k) = eos_smallt
          error(k)    = 0.1*eos_tol
       endif

    enddo
!!$    print*,'error  (NR1(Eos)):',error(vecBegin:vecEnd)

    ! Loop over the zones individually now
    do k = vecBegin, vecEnd
       do i = 2, eos_maxNewton
          if (error(k) .LT. eos_tol) goto 70

          call local_setEosTempsFromTempRow(k,k,ancMode,cMask)
!!$          if (.NOT. present(cMask)) then
!!$             tempRadRow(k) = tempRow(k)
!!$             tempIonRow(k) = tempRow(k)
!!$             tempEleRow(k) = tempRow(k)
!!$          else
!!$             tempRadRow(k) = eosdata(tempRad+k)
!!$             tempIonRow(k) = eosdata(tempIon+k)
!!$             tempEleRow(k) = eosdata(tempEle+k)
!!$             if (cMask(EOSCOMP_RAD).NE.0) tempRadRow(k) = tempRow(k)
!!$             if (cMask(EOSCOMP_ION).NE.0) tempIonRow(k) = tempRow(k)
!!$             if (cMask(EOSCOMP_ELE).NE.0) tempEleRow(k) = tempRow(k)
!!$          end if
!!$          select case (ancMode)
!!$          case(MODE_DENS_TEMP,MODE_DENS_TEMP_EQUI)
!!$             eosdata(temp+k) = tempRow(k)
!!$!          eosData(pres+k) = 0.0
!!$!             eosData(presIon+k) = 0.0!!!
!!$!             eosData(presEle+k) = 0.0
!!$!             eosData(presRad+k) = 0.0
!!$          case(MODE_DENS_TEMP_ALL)
!!$             eosdata(temp+k) = tempRow(k)
!!$             eosdata(tempIon+k) = tempIonRow(k)
!!$             eosdata(tempEle+k) = tempEleRow(k)
!!$             eosdata(tempRad+k) = tempRadRow(k)
!!$          case(MODE_DENS_TEMP_GATHER)
!!$             eosdata(tempIon+k) = tempIonRow(k)
!!$             eosdata(tempEle+k) = tempEleRow(k)
!!$             eosdata(tempRad+k) = tempRadRow(k)
!!$          case(MODE_DENS_TEMP_ION)
!!$             eosdata(tempIon+k) = tempIonRow(k)
!!$          case(MODE_DENS_TEMP_ELE)
!!$             eosdata(tempEle+k) = tempEleRow(k)
!!$          case(MODE_DENS_TEMP_RAD)
!!$             eosdata(tempRad+k) = tempRadRow(k)
!!$          end select
          ! do eos only over this single item
!!$          call eos_multiTypeByTemp(k,k,mask=eosMaskPlus,componentMask=cMask)
          call eos_multiTypeByTemp(ancMode,vecLen,eosData,k,k,massFrac,mask=eosMaskPlus,componentMask=cMask)
          call local_setTempRow(k,k,ancMode)
!!$          select case (ancMode)
!!$          case(MODE_DENS_TEMP_ION)
!!$             tempRow(k) = eosData(tempIon+k)
!!$          case(MODE_DENS_TEMP_ELE)
!!$             tempRow(k) = eosData(tempEle+k)
!!$          case(MODE_DENS_TEMP_RAD)
!!$             tempRow(k) = eosData(tempRad+k)
!!$          case default
!!$             eosData(presIon+vecBegin:presIon+vecEnd) = 0.0!!!
!!$             tempRow(k) = eosData(temp+k) !maybe not..
!!$          end select

          call local_setTots(k,k,mode)
!!$          select case (mode)
!!$          case(MODE_DENS_EI_ION,MODE_DENS_PRES_ION)
!!$             etotRow(k) = eosData(eintIon+k)
!!$             ptotRow(k) = eosData(presIon+k)
!!$          case(MODE_DENS_EI_ELE,MODE_DENS_PRES_ELE,MODE_DENS_ENTR_ELE)
!!$             etotRow(k) = eosData(eintEle+k)
!!$             ptotRow(k) = eosData(presEle+k)
!!$             stotRow(k) = eosData(entrEle+k)
!!$          case(MODE_DENS_EI_RAD,MODE_DENS_PRES_RAD)
!!$             etotRow(k) = eosData(eintRad+k)
!!$             ptotRow(k) = eosData(presRad+k)
!!$          case(MODE_DENS_EI)
!!$             etotRow(k) = eosData(eint+k)
!!$             ptotRow(k) = 0.0!!!
!!$          case default
!!$             etotRow(k) = eosData(eint+k)
!!$             ptotRow(k) = eosData(pres+k)
!!$          end select
          if (mode==MODE_DENS_EI)  ptotRow(k) = 0.0!!!  - why???

!!$          gamcRow(k) = eosData(gamc+k)
!!!$!!!          stotRow(k) = eosData(entr+k)
!!$          print*,'ptotRow(NR2(Eos)):',ptotRow(k),k,tempRow(k),etotRow(k),gamcRow(k)

          if (whatWanted == WANT_ENER) then
             if (eosData(det+k) .LE. 0.0) then
                print*,'Eos: negative dE/dT is bad!', eosData(det+k),tnew(k),tempRow(k),  ewantRow(k),etotRow(k)
             end if
             tnew(k) = tempRow(k) - (etotRow(k) - ewantRow(k)) / eosData(det+k)
          else if (whatWanted == WANT_ENTR) then
             tnew(k) = tempRow(k) - (stotRow(k) - swantRow(k)) / eosData(dst+k)
          else                    !pressure is wanted
             tnew(k) = tempRow(k) - (ptotRow(k) - pwantRow(k)) / eosData(dpt+k)
          end if
          ! if we jump out of the brackets, check if there is a solution
          if ( tnew(k) < bnd_lo(k) ) then
             told = tempRow(k)
             tempRow(k) = bnd_lo(k)
             call local_setEosTempsFromTempRow(k,k,ancMode,cMask)
             call eos_multiTypeByTemp(ancMode,vecLen,eosData,k,k,massFrac,mask=eosMaskPlus,componentMask=cMask)
             call local_setTempRow(k,k,ancMode)
             call local_setTots(k,k,mode)
             ! assume eint is a monotonically increasing function of T
             if ( xXwantRow(k) <= xXtotRow(k) ) then
                ! this will force exit from the loop via error=0
                tnew(k) = tempRow(k)
             else
                ! bisection
                tempRow(k) = told
                bnd_hi(k) = told
                tnew(k) = 0.5*(bnd_lo(k)+bnd_hi(k))
             endif
          else if ( tnew(k) .GE. bnd_hi(k) ) then
             told = tempRow(k)
             tempRow(k) = bnd_hi(k)
             call local_setEosTempsFromTempRow(k,k,ancMode,cMask)
             call eos_multiTypeByTemp(ancMode,vecLen,eosData,k,k,massFrac,mask=eosMaskPlus,componentMask=cMask)
             call local_setTempRow(k,k,ancMode)
             call local_setTots(k,k,mode)
             ! assume eint is a monotonically increasing function of T
             if ( xXwantRow(k) >= xXtotRow(k) ) then
                ! this will force exit from the loop via error=0
                tnew(k) = tempRow(k)
             else
                ! bisection
                tempRow(k) = told
                bnd_lo(k) = told
                tnew(k) = 0.5*(bnd_lo(k)+bnd_hi(k))
             endif
          else
             ! using derivative as an indicator of whether the old guess was
             ! an upper lower bound, update bracket range
             ! (assumes derivative is always well-behaved)
             if ( tnew(k) >= tempRow(k) ) then
                bnd_lo(k) = tempRow(k)
             else
                bnd_hi(k) = tempRow(k)
             endif
          endif

          ! Don't allow the temperature to change by more than an order of magnitude 
          ! in a single iteration
          if (tnew(k) .GT. eos_maxFactorUp*tempRow(k)) tnew(k) =  & 
               &              eos_maxFactorUp*tempRow(k)
          if (tnew(k) .LT. eos_maxFactorDown*tempRow(k)) tnew(k) =  & 
               &              eos_maxFactorDown*tempRow(k)

          ! Compute the error
          error(k) = abs((tnew(k) - tempRow(k)) / tempRow(k))
#ifdef DEBUG_EOS
          if (whatWanted == WANT_ENTR) then
             print*,'E2:',error(k),(stotRow(k) - swantRow(k)),tnew(k),tempRow(k) 
          else
          print*,'error  (NR2(Eos)):',error(k),tnew(k),tempRow(k)
          endif
#endif

          ! Store the new temperature
          tempRow(k) = tnew(k)

       ! DEV: Disable the following 'freezing' logic? - KW
          ! Check if we are freezing, if so set the temperature to eos_smallt, and adjust 
          ! the error so we don't wait for this one
          if (tempRow(k) .LT. eos_smallt) then
             tempRow(k) = eos_smallt
             error(k)    = .1*eos_tol
          endif

       end do  ! end of Newton iterations loop.  Failure drops below, success goes to 70

       ! Land here if too many iterations are needed -- failure

       print *, ' '
       print *, 'Newton-Raphson failed in multitemp, multitype Eos.F90'

       select case(mode)
       case(MODE_DENS_EI)
          print *, '(e and rho as input):'
       case(MODE_DENS_PRES)
          print *, '(pres and rho as input):'

       case(MODE_DENS_EI_ELE)
          print '(a)', "MODE_DENS_EI_ELE"
          print '(a)', "INPUTS: "
          print '(a,1pe24.15)', "mass density                      = ", denRow(k)
          print '(a,1pe24.15)', "electron specific internal energy = ", ewantRow(k)
          print '(/,a)', "CURRENT ITERATION:"
          print '(a,1pe24.15)', "electron specific internal energy = ", etotRow(k)
          print '(a,1pe24.15)', "electron temperature              = ", tempRow(k)

          print '(a)', "mass fractions:"
          do i = (k-1)*NSPECIES + 1, k*NSPECIES
             spec_num = i-(k-1)*NSPECIES
             call Simulation_mapIntToStr(spec_num+SPECIES_BEGIN-1,spec_str,MAPBLOCK_UNK)
             print '(a6,1pe24.15)', trim(spec_str), massfrac(i)
          end do

       case(MODE_DENS_ENTR_ELE)
          print '(a)', "MODE_DENS_ENTR_ELE"
          print '(a)', "INPUTS: "
          print '(a,1pe24.15)', "mass density                      = ", denRow(k)
          print '(a,1pe24.15)', "electron specific entropy         = ", swantRow(k)
          print '(/,a)', "CURRENT ITERATION:"
          print '(a,1pe24.15)', "electron specific entropy         = ", stotRow(k)
          print '(a,1pe24.15)', "electron specific internal energy = ", etotRow(k)
          print '(a,1pe24.15)', "electron temperature              = ", tempRow(k)

          print '(a)', "mass fractions:"
          do i = (k-1)*NSPECIES + 1, k*NSPECIES
             spec_num = i-(k-1)*NSPECIES
             call Simulation_mapIntToStr(spec_num+SPECIES_BEGIN-1,spec_str,MAPBLOCK_UNK)
             print '(a6,1pe24.15)', trim(spec_str), massfrac(i)
          end do

       case(MODE_DENS_EI_ION)
          print '(a)', "MODE_DENS_EI_ION"
          print '(a)', "INPUTS: "
          print '(a,1pe24.15)', "mass density                 = ", denRow(k)
          print '(a,1pe24.15)', "ion specific internal energy = ", ewantRow(k)
          print '(/,a)', "CURRENT ITERATION:"
          print '(a,1pe24.15)', "ion specific internal energy = ", etotRow(k)
          print '(a,1pe24.15)', "ion temperature              = ", tempRow(k)

          print '(a)', "mass fractions:"
          do i = (k-1)*NSPECIES + 1, k*NSPECIES
             spec_num = i-(k-1)*NSPECIES
             call Simulation_mapIntToStr(spec_num+SPECIES_BEGIN-1,spec_str,MAPBLOCK_UNK)
             print '(a6,1pe24.15)', trim(spec_str), massfrac(i)
          end do

       case default
          print *, 'Eos mode:', mode
       end select

       print *, ' '
       print *, 'too many iterations, reached', eos_maxNewton
       print *, ' '
       print *, ' k    = ', k,vecBegin,vecEnd
       print *, ' temp = ', tempRow(k)
       print *, ' dens = ', denRow(k)
       print *, ' pres = ', ptotRow(k)

       call Driver_abortFlash('[Eos] Error: too many iterations in Newton-Raphson')


       ! Land here if the Newton iteration converged
       !  jumps out of the iterations, but then continues to the next vector location

70     continue           
#ifdef DEBUG_EOS
       print*,'Eos nr done after',i,' iterations for',k
#endif
    end do

    ! Crank through the entire eos one last time

    call local_setEosTempsFromTempRow(ilo,ihi,ancMode,cMask)
!!$    if (.NOT. present(cMask)) then
!!$       tempRadRow(ilo:ihi) = tempRow(ilo:ihi)
!!$       tempIonRow(ilo:ihi) = tempRow(ilo:ihi)
!!$       tempEleRow(ilo:ihi) = tempRow(ilo:ihi)
!!$    else
!!$       tempRadRow(ilo:ihi) = eosdata(tempRad+ilo:tempRad+ihi)
!!$       tempIonRow(ilo:ihi) = eosdata(tempIon+ilo:tempIon+ihi)
!!$       tempEleRow(ilo:ihi) = eosdata(tempEle+ilo:tempEle+ihi)
!!$       if (cMask(EOSCOMP_RAD).NE.0) tempRadRow(ilo:ihi) = tempRow(ilo:ihi)
!!$       if (cMask(EOSCOMP_ION).NE.0) tempIonRow(ilo:ihi) = tempRow(ilo:ihi)
!!$       if (cMask(EOSCOMP_ELE).NE.0) tempEleRow(ilo:ihi) = tempRow(ilo:ihi)
!!$    end if
!!$    select case (ancMode)
!!$    case(MODE_DENS_TEMP,MODE_DENS_TEMP_EQUI) 
!!$       eosdata(temp+ilo:temp+ihi)       = tempRow(ilo:ihi)
!!$!          eosData(pres+k) = 0.0
!!$!             eosData(presIon+k) = 0.0
!!$!             eosData(presEle+k) = 0.0
!!$!             eosData(presRad+k) = 0.0
!!$!       eosData(presIon+vecBegin:presIon+vecEnd) = 123.45
!!$    case(MODE_DENS_TEMP_ALL)
!!$       eosdata(temp+ilo:temp+ihi)       = tempRow(ilo:ihi)
!!$       eosdata(tempIon+ilo:tempIon+ihi) = tempIonRow(ilo:ihi)
!!$       eosdata(tempEle+ilo:tempEle+ihi) = tempEleRow(ilo:ihi)
!!$       eosdata(tempRad+ilo:tempRad+ihi) = tempRadRow(ilo:ihi)
!!$    case(MODE_DENS_TEMP_GATHER)
!!$       eosdata(tempIon+ilo:tempIon+ihi) = tempIonRow(ilo:ihi)
!!$       eosdata(tempEle+ilo:tempEle+ihi) = tempEleRow(ilo:ihi)
!!$       eosdata(tempRad+ilo:tempRad+ihi) = tempRadRow(ilo:ihi)
!!$    case(MODE_DENS_TEMP_ION)
!!$       eosdata(tempIon+ilo:tempIon+ihi) = tempIonRow(ilo:ihi)
!!$    case(MODE_DENS_TEMP_ELE)
!!$       eosdata(tempEle+ilo:tempEle+ihi) = tempEleRow(ilo:ihi)
!!$    case(MODE_DENS_TEMP_RAD)
!!$       eosdata(tempRad+ilo:tempRad+ihi) = tempRadRow(ilo:ihi)
!!$    end select

!!$    call eos_multiTypeByTemp(vecBegin,vecEnd,mask=eosMaskPlus,componentMask=cMask)
    call eos_multiTypeByTemp(ancMode,vecLen,eosData,vecBegin,vecEnd,massFrac,mask=eosMaskPlus,componentMask=cMask)
    tempRow(vecBegin:vecEnd) = eosData(temp+vecBegin:temp+vecEnd) !maybe not..
    call local_setTots(vecBegin,vecEnd,mode)
!!$    select case (mode)
!!$    case(MODE_DENS_EI_ION,MODE_DENS_PRES_ION)
!!$       etotRow(vecBegin:vecEnd) = eosData(eintIon+vecBegin:eintIon+vecEnd)
!!$       ptotRow(vecBegin:vecEnd) = eosData(presIon+vecBegin:presIon+vecEnd)
!!$    case(MODE_DENS_EI_ELE,MODE_DENS_PRES_ELE)
!!$       etotRow(vecBegin:vecEnd) = eosData(eintEle+vecBegin:eintEle+vecEnd)
!!$       ptotRow(vecBegin:vecEnd) = eosData(presEle+vecBegin:presEle+vecEnd)
!!$    case(MODE_DENS_EI_RAD,MODE_DENS_PRES_RAD)
!!$       etotRow(vecBegin:vecEnd) = eosData(eintRad+vecBegin:eintRad+vecEnd)
!!$       ptotRow(vecBegin:vecEnd) = eosData(presRad+vecBegin:presRad+vecEnd)
!!$    case(MODE_DENS_EI)
!!$       etotRow(vecBegin:vecEnd) = eosData(eint+vecBegin:eint+vecEnd)
!!$       ptotRow(vecBegin:vecEnd) = eosData(pres+vecBegin:pres+vecEnd)
!!$    case default
!!$       etotRow(vecBegin:vecEnd) = eosData(eint+vecBegin:eint+vecEnd)
!!$       ptotRow(vecBegin:vecEnd) = eosData(pres+vecBegin:pres+vecEnd)
!!$    end select

!!$    etotRow(vecBegin:vecEnd) = eosData(eint+vecBegin:eint+vecEnd)
!!$    ptotRow(vecBegin:vecEnd) = eosData(pres+vecBegin:pres+vecEnd)
!!$    gamcRow(vecBegin:vecEnd) = eosData(gamc+vecBegin:gamc+vecEnd)
!!!$!!!    stotRow(vecBegin:vecEnd) = eosData(entr+vecBegin:entr+vecEnd)
!!$    print*,'ptotRow(NR3(Eos)):',ptotRow(vecBegin:vecEnd)

    ! Fill the FLASH arrays with the results.  
    !  In MODE_DENS_EI, we should be generating temperature and pressure (plus gamma and entropy)

    select case(mode)
    case(MODE_DENS_EI)
       eosData(temp+ilo:temp+ihi)=tempRow(ilo:ihi)
!       eosData(pres+ilo:pres+ihi)=ptotRow(ilo:ihi) ! moved up to callers - KW 2010-12-13
       !  Update the energy to be the true energy, instead of the energy we were trying to meet
       !  ConstantInput LBR and KW believe this is WRONG -- the input arrays should not be changed
       if (forceConstInp)  then
          eosData(eint+ilo:eint+ihi) = esaveRow(ilo:ihi)
       else
          eosData(eint+ilo:eint+ihi) = etotRow(ilo:ihi)
       end if
    case(MODE_DENS_EI_ION)
       if (forceConstInp)  then
          eosData(eintIon+ilo:eintIon+ihi) = esaveRow(ilo:ihi)
       else
          eosData(eintIon+ilo:eintIon+ihi) = etotRow(ilo:ihi)
       end if
    case(MODE_DENS_EI_ELE)
       if (forceConstInp)  then
          eosData(eintEle+ilo:eintEle+ihi) = esaveRow(ilo:ihi)
       else
          eosData(eintEle+ilo:eintEle+ihi) = etotRow(ilo:ihi)
       end if
    case(MODE_DENS_EI_RAD)
       if (forceConstInp)  then
          eosData(eintRad+ilo:eintRad+ihi) = esaveRow(ilo:ihi)
       else
          eosData(eintRad+ilo:eintRad+ihi) = etotRow(ilo:ihi)
       end if
    case(MODE_DENS_PRES)
       eosData(temp+ilo:temp+ihi)=tempRow(ilo:ihi)
       ! Update the pressure to be the equilibrium pressure, instead of the pressure we were trying to meet
       !  ConstantInput LBR and KW believe this is wrong.  See notes at the top of the routine
       if (forceConstInp) then
          eosData(pres+ilo:pres+ihi) = psaveRow(ilo:ihi)
       else
          eosData(pres+ilo:pres+ihi) = ptotRow(ilo:ihi)
       end if
       eosData(eint+ilo:eint+ihi)=etotRow(ilo:ihi)
    case(MODE_DENS_PRES_ION)
       if (forceConstInp)  then
          eosData(presIon+ilo:presIon+ihi) = psaveRow(ilo:ihi)
       else
          eosData(presIon+ilo:presIon+ihi) = ptotRow(ilo:ihi)
       end if
    case(MODE_DENS_PRES_ELE)
       if (forceConstInp)  then
          eosData(presEle+ilo:presEle+ihi) = psaveRow(ilo:ihi)
       else
          eosData(presEle+ilo:presEle+ihi) = ptotRow(ilo:ihi)
       end if
    case(MODE_DENS_PRES_RAD)
       if (forceConstInp)  then
          eosData(presRad+ilo:presRad+ihi) = psaveRow(ilo:ihi)
       else
          eosData(presRad+ilo:presRad+ihi) = ptotRow(ilo:ihi)
       end if
    end select

!    eosData(gamc+ilo:gamc+ihi)=gamcRow(ilo:ihi) ! moved up to callers - KW 2010-12-13
!    eosData(entr+ilo:entr+ihi)=stotRow(ilo:ihi) ! moved up to callers - KW 2010-12-13
#ifdef DEBUG_DEC2010
    print*,'NR+I',eCompRow(EOSCOMP_ION,1:vecLen)
!    print*,'NR+e',eCompRow(EOSCOMP_ELE,1:vecLen)
!    print*,'NR+R',eCompRow(EOSCOMP_RAD,1:vecLen)
#endif

  contains
    subroutine local_setEosTempsFromTempRow(a,b,ancMode,cMask)
      integer,intent(IN) :: a,b,ancMode
      integer,optional, dimension(EOSCOMP_NUM_COMPONENTS),INTENT(in)::cMask
      if (.NOT. present(cMask)) then
         tempRadRow(a:b) = tempRow(a:b)
         tempIonRow(a:b) = tempRow(a:b)
         tempEleRow(a:b) = tempRow(a:b)
      else
         tempRadRow(a:b) = eosdata(tempRad+a:tempRad+b)
         tempIonRow(a:b) = eosdata(tempIon+a:tempIon+b)
         tempEleRow(a:b) = eosdata(tempEle+a:tempEle+b)
         if (cMask(EOSCOMP_RAD).NE.0) tempRadRow(a:b) = tempRow(a:b)
         if (cMask(EOSCOMP_ION).NE.0) tempIonRow(a:b) = tempRow(a:b)
         if (cMask(EOSCOMP_ELE).NE.0) tempEleRow(a:b) = tempRow(a:b)
      end if
      select case (ancMode)
      case(MODE_DENS_TEMP,MODE_DENS_TEMP_EQUI)
         eosdata(temp+a:temp+b)       = tempRow(a:b)
      case(MODE_DENS_TEMP_ALL)
         eosdata(temp+a:temp+b)       = tempRow(a:b)
         eosdata(tempIon+a:tempIon+b) = tempIonRow(a:b)
         eosdata(tempEle+a:tempEle+b) = tempEleRow(a:b)
         eosdata(tempRad+a:tempRad+b) = tempRadRow(a:b)
      case(MODE_DENS_TEMP_GATHER)
         eosdata(tempIon+a:tempIon+b) = tempIonRow(a:b)
         eosdata(tempEle+a:tempEle+b) = tempEleRow(a:b)
         eosdata(tempRad+a:tempRad+b) = tempRadRow(a:b)
      case(MODE_DENS_TEMP_ION)
         eosdata(tempIon+a:tempIon+b) = tempIonRow(a:b)
      case(MODE_DENS_TEMP_ELE)
         eosdata(tempEle+a:tempEle+b) = tempEleRow(a:b)
      case(MODE_DENS_TEMP_RAD)
         eosdata(tempRad+a:tempRad+b) = tempRadRow(a:b)
      end select
    end subroutine local_setEosTempsFromTempRow

    subroutine local_setTempRow(a,b,ancMode)
      integer,intent(IN) :: a,b,ancMode
      select case (ancMode)
      case(MODE_DENS_TEMP_ION)
         tempRow(a:b) = eosData(tempIon+a:tempIon+b)
      case(MODE_DENS_TEMP_ELE)
         tempRow(a:b) = eosData(tempEle+a:tempEle+b)
      case(MODE_DENS_TEMP_RAD)
         tempRow(a:b) = eosData(tempRad+a:tempRad+b)
      case default
         tempRow(a:b) = eosData(temp+a:temp+b)
      end select
    end subroutine local_setTempRow
    subroutine local_setTots(a,b,mode)
      integer,intent(IN) :: a,b,mode
      select case (mode)
      case(MODE_DENS_EI_ION,MODE_DENS_PRES_ION)
         etotRow(a:b) = eosData(eintIon+a:eintIon+b)
         ptotRow(a:b) = eosData(presIon+a:presIon+b)
      case(MODE_DENS_EI_ELE,MODE_DENS_PRES_ELE)
         etotRow(a:b) = eosData(eintEle+a:eintEle+b)
         ptotRow(a:b) = eosData(presEle+a:presEle+b)
      case(MODE_DENS_ENTR_ELE)
         etotRow(a:b) = eosData(eintEle+a:eintEle+b)
         ptotRow(a:b) = eosData(presEle+a:presEle+b)
         stotRow(a:b) = eosData(entrEle+a:entrEle+b)
      case(MODE_DENS_EI_RAD,MODE_DENS_PRES_RAD)
         etotRow(a:b) = eosData(eintRad+a:eintRad+b)
         ptotRow(a:b) = eosData(presRad+a:presRad+b)
      case default
         etotRow(a:b) = eosData(eint+a:eint+b)
         ptotRow(a:b) = eosData(pres+a:pres+b)
      end select

    end subroutine local_setTots
  end subroutine eos_newtonRaphson

