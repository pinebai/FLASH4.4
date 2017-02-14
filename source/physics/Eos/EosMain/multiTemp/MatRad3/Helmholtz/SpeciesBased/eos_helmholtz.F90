!!****if* source/physics/Eos/EosMain/multiTemp/MatRad3/Helmholtz/SpeciesBased/eos_helmholtz
!!
!! NAME
!!
!! eos_helmholtz
!!
!! SYNOPSIS
!!
!!  call eos_helmholtz(integer(IN) :: mode,
!!                     integer(IN) :: vecLen,
!!                     real(INOUT) :: eosData(vecLen*EOS_NUM),
!!           optional, real(IN)    :: massFrac(vecLen*NSPECIES),
!!     optional,target,logical(IN) :: mask(EOS_VARS+1:EOS_NUM)    )
!!
!! DESCRIPTION
!!
!!   Driver for the Helmholtz and Nadyozhin equations of state.
!!   See the NOTES section for important information about this implementation.
!!
!!  This routine applies the equation of state to thermodynamic 
!!  quantities at one or more grid cells.  The number of cells is 
!!  determined by the argument veclen.  Data is packaged for this 
!!  routine in the 1d array, eosData.  The data in eosData is organized: 
!!  1:vecLen points contain the first variable, vecLen+1:2*vecLen points 
!!  contain the second variable, and so on. The number and order of
!!  variables in the array is determined by the constants defined in Eos.h.
!!  
!!  The routine takes different quantities as givens depending on the
!!  value of the mode variable: if mode=MODE_DENS_TEMP, density and
!!  temperature are taken as given, and pressure and energy are generated
!!  as output; if mode=MODE_DENS_EI, density and internal energy are taken as
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
!!  ARGUMENTS
!! 
!!  mode :    Selects the mode of operation of the Eos unit.
!!             The valid values are MODE_DENS_EI, MODE_DENS_PRES and  
!!             MODE_DENS_TEMP as decribed above.
!!
!!  vecLen   : number of points (cells) for which the eosData array is sized.
!!             This is the
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
!!             call without a mask actual acgument, or set the mask equal to .false.
!!
!!             The Helmholtz EOS kernel calculation ignores the mask setting and calculates
!!             all derivatives, whether needed or not.  This routine does not return
!!             derivatives if the mask is requested, but the calculation is not speeded up
!!             significantly by setting the mask.
!!
!!
!! PARAMETERS
!!
!!  eos_tol    Controls the accuracy of the Newton-Raphson iterations for MODE_DENS_EI and 
!!             MODE_DENS_PRES.
!!
!! EXAMPLE
!!
!! --- A single-point at a time example, does not calculate derivatives (based on Cellular Simulation)---
!!
!!  #include "constants.h"   ! for MODE_DENS_TEMP
!!  #include "Flash.h"       ! for NSPECIES
!!  #include "Eos.h"         ! for EOS_VAR order
!!
!!  real  :: temp_zone, rho_zone, ptot, eint, entr, gamma
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
!!           entr = eosData(EOS_ENTR)
!!           gamma = eosData(EOS_GAMC)
!!           
!!           call Grid_putPointData(blockId,CENTER,TEMP_VAR,EXTERIOR,iPosition,temp_zone)
!!           call Grid_putPointData(blockId,CENTER,DENS_VAR,EXTERIOR,iPosition,rho_zone)
!!           call Grid_putPointData(blockId,CENTER,PRES_VAR,EXTERIOR,iPosition,ptot)
!!           call Grid_putPointData(blockId,CENTER,EINT_VAR,EXTERIOR,iPosition,eint)
!!               if you want ENER_VAR, calculate it from eint and kinetic energy
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
!!  The Helmholtz equation of state calculations are iterative for any mode other
!!  than MODE_DENS_TEMP.  Therefore, the intial estimates for temperature and density
!!  must be pretty good upon entering Eos with any other MODE_....or the calculations will
!!  not converge.
!!
!!  This algorithm uses a data table helm_table.dat which contains the coefficients for
!!  one of the interpolating algorithms.  Upon first entry to the Eos, a binary version of this
!!  table (helm_table.bdat) is created for speed of access.  This binary file should NOT be
!!  carried across machine architectures or even compilers.
!!
!!  When USE_EOS_YE is defined, the implementation in
!!  physics/Eos/EosMain/multiTemp/Helmholtz/Ye is supposed to be used rather than this
!!  species-based one.
!!
!!  When operating in MODE_DENS_EI, the INPUT energy is updated.  This change of an input parameter
!!     can be overridden by setting the runtime parameter eos_forceConstantInput to true.
!!     Noted below, see comments prefaced with ConstantInput.
!!  Similarly, when operating in MODE_DENS_PRES, the INPUT pressure is updated.  Physicists need
!!     to be aware of this.  Similarly can be overridden with the runtime parameter/
!!
!!  The accuracy can be adjusted with the parameter eos_tol.
!!
!! SEE ALSO
!! 
!!  Eos.h    defines the variables used.
!!  Eos_wrapped  sets up the data structure.
!!
!!
!!*** 


subroutine eos_helmholtz(mode,vecLen,eosData,massFrac,mask)

  use Driver_interface, ONLY : Driver_abortFlash
  use Multispecies_interface, ONLY : Multispecies_getSumInv, &
    Multispecies_getSumFrac
  use Logfile_interface, ONLY:  Logfile_stampMessage
  use Hydro_interface, ONLY : Hydro_recalibrateEints

  use eos_helmInterface, ONLY : eos_helm

  use eos_helmData, ONLY: eos_tol, eos_maxNewton
  use Eos_data, ONLY : eos_smallT, eos_meshMe, eos_singleSpeciesA, eos_singleSpeciesZ,&
       eos_combinedTempRule,eos_forceConstantInput, eos_largeT
  use eos_helmConstData, ONLY: eos_ao3 !, eos_kerg, eos_avo, eos_kergavo, eos_h,eos_hbar
  use eos_vecData, ONLY:  tempRow, denRow, etotRow, abarRow, zbarRow, &
       gamcRow, ptotRow, deaRow, dezRow, stotRow, dsdRow, dstRow, &
       detRow, dptRow, dpdRow, dedRow, pelRow, neRow, etaRow, cvRow, cpRow, &
       cvIonRow, cvEleRow, &
       tempRadRow,tempIonRow,tempEleRow,eCompRow,pCompRow

  use Simulation_interface, ONLY: Simulation_mapIntToStr

  implicit none

#include "constants.h"
#include "Flash.h"
#include "Eos.h"
#include "Eos_components.h"
#ifdef FLASH_MULTISPECIES
#include "Multispecies.h"
#endif

  !     Arguments
  integer, INTENT(in) :: mode, vecLen
  real, INTENT(inout), dimension(vecLen*EOS_NUM) :: eosData
  real, optional,INTENT(in), dimension(vecLen*NSPECIES) :: massFrac
  ! must correspond to dimensions of Eos_wrapped
  logical, optional, INTENT(in),target,dimension(EOS_VARS+1:EOS_NUM) :: mask

  ! This is the variable that is used internally -- set to false unless mask comes in
  logical,save,target, dimension(EOS_VARS+1:EOS_NUM) :: maskInternal
  data maskInternal / EOS_DERIVS * .FALSE. /
  logical,pointer, dimension(:) :: maskPtr
  !! MANOS B
!!$  logical, dimension(EOS_VARS+1:EOS_NUM) :: maskInternal
  !! MANOS E

  integer :: componentMask(EOSCOMP_NUM_COMPONENTS)
  real,parameter :: eos_gammaRad = (4./3.)

  integer :: i, k
  integer :: vecBegin,vecEnd
  integer :: pres, temp, dens, gamc, eint
  integer :: tempIon,tempEle,tempRad
  integer :: eintIon,eintEle,eintRad
  integer :: presIon,presEle,presRad
  integer :: entrEle,entrRad
  integer :: abar, zbar
  integer :: entr, dst, dsd
  integer :: dpt, dpd, det, ded, dea, dez, pel, ne, eta, c_v, c_p
  real    :: abarInv, zbarFrac
  integer :: ilo,ihi, rowLen
  logical :: arrayBoundHackMode

  ! declare some local storage for the results of the Newton iteration
  real,dimension(vecLen)::  ewantRow, tnew, error,pwantRow
  target :: ewantRow,pwantRow
  !  local storage for forcedConstantInput -- could be allocatable, but might be so slow
  !  that it's not worth the small storage save.
  real,dimension(vecLen)::  psaveRow, esaveRow

  integer :: spec_num
  character(len=MAX_STRING_LENGTH) :: spec_str

  logical  :: firstCall = .true.

  !      Fill the pipe with the initial temperature, density, and composition.
  !      The composition is parametrized by abar and zbar, which are computed
  !      from the mass fractions xn().

  !  if you are working with electron abundance mass scalars, then you don't
  !  necessarily have to have mass fractions.
  if(.not.present(massFrac)) then
     call Driver_abortFlash("[Eos] Helmholtz needs mass fractions")
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

  vecBegin = 1
  vecEnd = vecLen

  ilo = vecBegin
  ihi = vecEnd
  rowLen = ihi - ilo + 1

  ! These integers are indexes into the lowest location in UNK that contain the appropriate variable
  pres = (EOS_PRES-1)*vecLen
  dens = (EOS_DENS-1)*vecLen
  temp = (EOS_TEMP-1)*vecLen
  eint = (EOS_EINT-1)*vecLen   ! in flash2 eos_helm, this is etot
  gamc = (EOS_GAMC-1)*vecLen   ! in flash2 eos_helm, this is gamc
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

  !! For allocatable arrays, set them up now.
#ifndef FIXEDBLOCKSIZE
  call eos_vecAlloc(vecLen)
#endif

  do k = 1, vecLen

     tempRow(k)    = eosData(temp+k)
     denRow(k)     = eosData(dens+k)
     ewantRow(k)   = eosData(eint+k)   ! store desired internal energy for mode=2 case
     pwantRow(k)   = eosData(pres+k)   ! store desired pressure for mode=3 case

     ! Note in Eos.F90, we assume the user knows what he's doing.  Eos_wrapped does not.

#ifdef FLASH_MULTISPECIES
     !Calculate the inverse in a way that allows for zero mass fractions
     call Multispecies_getSumInv(A, abarInv,massFrac((k-1)*NSPECIES+1:k*NSPECIES))
     abarRow(k) = 1.e0 / abarInv

     call Multispecies_getSumFrac(Z, zbarFrac, massFrac((k-1)*NSPECIES+1:k*NSPECIES))
     zbarRow(k) = abarRow(k) * zbarFrac
#else
     ! No multispecies defined, use default values (same as Gamma formulation)
     abarRow(k) = eos_singleSpeciesA
     zbarRow(k) = eos_singleSpeciesZ
#endif

     if (mode==MODE_DENS_EI_RECAL_GATHER) then
        call Hydro_recalibrateEints(eosData(eint+k), &
             eosData(eintIon+k),eosData(eintEle+k),eosData(eintRad+k))
     end if
  enddo

  !! Save the input parameters if eos_forceConstantInput is defined
  if (eos_forceConstantInput) then
     esaveRow = ewantRow
     psaveRow = pwantRow
  end if


  eosData(abar+1:abar+vecLen) = abarRow(1:vecLen) 
  eosData(zbar+1:zbar+vecLen) = zbarRow(1:vecLen)

  !==============================================================================
  !      MODE_DENS_TEMP  temperature and density given

  !      Crank the EOS on the pipes filled above, then fill the FLASH arrays
  !      with the thermodynamic quantities returned by the EOS.

  select case (mode)
  case (MODE_DENS_TEMP)

     tempRadRow(1:vecLen) = tempRow(1:vecLen)
     tempIonRow(1:vecLen) = tempRow(1:vecLen)
     tempEleRow(1:vecLen) = tempRow(1:vecLen)

     call eos_helm(1,vecLen,maskPtr)

     eosData(pres+1:pres+vecLen)=ptotRow(1:vecLen)
     eosData(eint+1:eint+vecLen)=etotRow(1:vecLen)
     eosData(gamc+1:gamc+vecLen)=gamcRow(1:vecLen)
     eosData(entr+1:entr+vecLen)=stotRow(1:vecLen)

  case(MODE_DENS_TEMP_ION,MODE_DENS_TEMP_ELE,MODE_DENS_TEMP_RAD)

     componentMask(:) = 0
     select case (mode)
     case(MODE_DENS_TEMP_RAD)
        componentMask(EOSCOMP_RAD) = 1
     case(MODE_DENS_TEMP_ELE)
        componentMask(EOSCOMP_ELE) = 1
     case(MODE_DENS_TEMP_ION)
        componentMask(EOSCOMP_ION) = 1
     end select
     tempRadRow(1:vecLen) = eosData(tempRad+1:tempRad+vecLen)
     tempIonRow(1:vecLen) = eosData(tempIon+1:tempIon+vecLen)
     tempEleRow(1:vecLen) = eosData(tempEle+1:tempEle+vecLen)

     call eos_helm(1,vecLen,maskPtr,componentMask)

     if(maskPtr(EOS_PRESION)) eosData(presIon+1:presIon+vecLen)=0.0
     if(maskPtr(EOS_PRESELE)) eosData(presEle+1:presEle+vecLen)=pCompRow(EOSCOMP_ELE,1:vecLen)+pCompRow(EOSCOMP_ION,1:vecLen)
     if(maskPtr(EOS_PRESRAD)) eosData(presRad+1:presRad+vecLen)=pCompRow(EOSCOMP_RAD,1:vecLen)
     if(maskPtr(EOS_EINTION)) eosData(eintIon+1:eintIon+vecLen)=0.0
     if(maskPtr(EOS_EINTELE)) eosData(eintEle+1:eintEle+vecLen)=eCompRow(EOSCOMP_ELE,1:vecLen)+eCompRow(EOSCOMP_ION,1:vecLen)
     if(maskPtr(EOS_EINTRAD)) eosData(eintRad+1:eintRad+vecLen)=eCompRow(EOSCOMP_RAD,1:vecLen)
     eosData(pres+1:pres+vecLen)=ptotRow(1:vecLen)
     eosData(eint+1:eint+vecLen)=etotRow(1:vecLen)
     eosData(gamc+1:gamc+vecLen)=gamcRow(1:vecLen)
     eosData(entr+1:entr+vecLen)=stotRow(1:vecLen)

  case(MODE_DENS_TEMP_COMP,MODE_DENS_TEMP_GATHER)

     tempRow(1:vecLen)    = eosData(temp+1:temp+vecLen)
     tempRadRow(1:vecLen) = eosData(tempRad+1:tempRad+vecLen)
     tempIonRow(1:vecLen) = eosData(tempEle+1:tempEle+vecLen)
     tempEleRow(1:vecLen) = eosData(tempEle+1:tempEle+vecLen)

     call eos_helm(1,vecLen,maskPtr)

     eosData(tempIon+1:tempIon+vecLen)= 0.0
     eosData(pres+1:pres+vecLen)=ptotRow(1:vecLen)
     if(maskPtr(EOS_PRESION)) eosData(presIon+1:presIon+vecLen)=0.0
     if(maskPtr(EOS_PRESELE)) eosData(presEle+1:presEle+vecLen)=pCompRow(EOSCOMP_ELE,1:vecLen)+pCompRow(EOSCOMP_ION,1:vecLen)
     if(maskPtr(EOS_PRESRAD)) eosData(presRad+1:presRad+vecLen)=pCompRow(EOSCOMP_RAD,1:vecLen)
     eosData(eint+1:eint+vecLen)=etotRow(1:vecLen)
     if(maskPtr(EOS_EINTION)) eosData(eintIon+1:eintIon+vecLen)=0.0
     if(maskPtr(EOS_EINTELE)) eosData(eintEle+1:eintEle+vecLen)=eCompRow(EOSCOMP_ELE,1:vecLen)+eCompRow(EOSCOMP_ION,1:vecLen)
     if(maskPtr(EOS_EINTRAD)) eosData(eintRad+1:eintRad+vecLen)=eCompRow(EOSCOMP_RAD,1:vecLen)
!!$     if(maskPtr(EOS_EINTELE)) print*,'eintEle:',eosData(eintEle+1:eintEle+vecLen)
!!$     if(maskPtr(EOS_EINTRAD)) print*,'eintRad:',eosData(eintRad+1:eintRad+vecLen)
     eosData(gamc+1:gamc+vecLen)=gamcRow(1:vecLen)
     eosData(entr+1:entr+vecLen)=stotRow(1:vecLen)

     if(mode == MODE_DENS_TEMP_GATHER)then
        if(eos_combinedTempRule==0) then
           call eos_newtonRaphson(vecLen, MODE_DENS_EI, .TRUE., maskPtr,cMask=(/1,1,1/))
        endif
        call setCombinedTemp()
     endif

  case (MODE_DENS_TEMP_EQUI)

     tempRadRow(1:vecLen) = tempRow(1:vecLen)
     tempIonRow(1:vecLen) = tempRow(1:vecLen)
     tempEleRow(1:vecLen) = tempRow(1:vecLen)

     call eos_helm(1,vecLen,maskPtr)

     eosData(tempIon+1:tempIon+vecLen)= 0.0
     eosData(tempEle+1:tempEle+vecLen)= tempRow(1:vecLen)
     eosData(tempRad+1:tempRad+vecLen)= tempRow(1:vecLen)
     eosData(pres+1:pres+vecLen)      = ptotRow(1:vecLen)
     if(maskPtr(EOS_PRESION)) eosData(presIon+1:presIon+vecLen)= 0.0
     if(maskPtr(EOS_PRESELE)) eosData(presEle+1:presEle+vecLen)= pCompRow(EOSCOMP_ELE,1:vecLen)+pCompRow(EOSCOMP_ION,1:vecLen)
     if(maskPtr(EOS_PRESRAD)) eosData(presRad+1:presRad+vecLen)= pCompRow(EOSCOMP_RAD,1:vecLen)
     eosData(eint+1:eint+vecLen)      = etotRow(1:vecLen)
     if(maskPtr(EOS_EINTION)) eosData(eintIon+1:eintIon+vecLen)=0.0
     if(maskPtr(EOS_EINTELE)) eosData(eintEle+1:eintEle+vecLen)=eCompRow(EOSCOMP_ELE,1:vecLen)+eCompRow(EOSCOMP_ION,1:vecLen)
     if(maskPtr(EOS_EINTRAD)) eosData(eintRad+1:eintRad+vecLen)=eCompRow(EOSCOMP_RAD,1:vecLen)
     eosData(gamc+1:gamc+vecLen)      = gamcRow(1:vecLen)
     eosData(entr+1:entr+vecLen)      = stotRow(1:vecLen)

      if (.NOT. all(tempRow(1:vecLen) > 0)) &
           print*,'eos_nr WARN _EQUI:',tempRow(1:vecLen)

     !==============================================================================
     !      MODE_DENS_EI  internal energy and density given

  case (MODE_DENS_EI, MODE_DENS_EI_SCATTER)

     call eos_newtonRaphson(vecLen, MODE_DENS_EI, .TRUE., maskPtr,cMask=(/1,1,1/))

     if (mode == MODE_DENS_EI_SCATTER) then
        eosData(tempIon+1:tempIon+vecLen)=0.0
        eosData(tempEle+1:tempEle+vecLen)=tempRow(1:vecLen)
        eosData(tempRad+1:tempRad+vecLen)=tempRow(1:vecLen)
        if(maskPtr(EOS_PRESION)) eosData(presIon+1:presIon+vecLen)= 0.0
        if(maskPtr(EOS_PRESELE)) eosData(presEle+1:presEle+vecLen)= pCompRow(EOSCOMP_ELE,1:vecLen)+pCompRow(EOSCOMP_ION,1:vecLen)
        if(maskPtr(EOS_PRESRAD)) eosData(presRad+1:presRad+vecLen)= pCompRow(EOSCOMP_RAD,1:vecLen)
        if(maskPtr(EOS_EINTION)) eosData(eintIon+1:eintIon+vecLen)=0.0
        if(maskPtr(EOS_EINTELE)) eosData(eintEle+1:eintEle+vecLen)=eCompRow(EOSCOMP_ELE,1:vecLen)+eCompRow(EOSCOMP_ION,1:vecLen)
        if(maskPtr(EOS_EINTRAD)) eosData(eintRad+1:eintRad+vecLen)=eCompRow(EOSCOMP_RAD,1:vecLen)
     end if


  case(MODE_DENS_EI_ION,MODE_DENS_EI_ELE,MODE_DENS_EI_RAD)

     componentMask(:) = 0
     tempRadRow(1:vecLen) = max(eos_smallT,eosData(tempRad+1:tempRad+vecLen))
     tempIonRow(1:vecLen) = eosData(tempIon+1:tempIon+vecLen)
     tempEleRow(1:vecLen) = eosData(tempEle+1:tempEle+vecLen)
     select case (mode)
     case(MODE_DENS_EI_RAD)
        componentMask(EOSCOMP_RAD) = 1
        ewantRow(1:vecLen) = eosData(eintRad+1:eintRad+vecLen)   ! store desired radiation energy
        call eos_newtonRaphson(vecLen, mode, .TRUE., maskPtr,componentMask)
        eosData(tempRad+1:tempRad+vecLen)=tempRow(1:vecLen)
        if(maskPtr(EOS_PRESRAD)) eosData(presRad+1:presRad+vecLen)=pCompRow(EOSCOMP_RAD,1:vecLen)
     case(MODE_DENS_EI_ELE)
        componentMask(EOSCOMP_ELE) = 1
        ewantRow(1:vecLen) = eosData(eintEle+1:eintEle+vecLen)   ! store desired electron energy
        call eos_newtonRaphson(vecLen, mode, .TRUE., maskPtr,componentMask)
        eosData(tempEle+1:tempEle+vecLen)=tempRow(1:vecLen)
        if(maskPtr(EOS_PRESELE)) eosData(presEle+1:presEle+vecLen)=pCompRow(EOSCOMP_ELE,1:vecLen)
     case(MODE_DENS_EI_ION)
        componentMask(EOSCOMP_ION) = 1
        ewantRow(1:vecLen) = eosData(eintIon+1:eintIon+vecLen)   ! store desired ion energy
        call eos_newtonRaphson(vecLen, mode, .TRUE., maskPtr,componentMask)
        eosData(tempIon+1:tempIon+vecLen)=tempRow(1:vecLen)
        if(maskPtr(EOS_PRESION)) eosData(presIon+1:presIon+vecLen)=pCompRow(EOSCOMP_ION,1:vecLen)
     end select


  case (MODE_DENS_EI_COMP, MODE_DENS_EI_GATHER, MODE_DENS_EI_RECAL_GATHER, MODE_DENS_EI_ALL)

     eosData(presIon+ilo:presIon+ihi) = 0.0
     eosData(tempIon+ilo:tempIon+ihi) = 0.0

     tempRadRow(1:vecLen) = max(eos_smallT,eosData(tempRad+1:tempRad+vecLen))
     tempIonRow(1:vecLen) = eosData(tempEle+1:tempEle+vecLen)
     tempEleRow(1:vecLen) = eosData(tempEle+1:tempEle+vecLen)

     componentMask(:) = (/0,0,0/)
     componentMask(EOSCOMP_RAD) = 1
     tempRow(1:vecLen) = tempRadRow(1:vecLen)
     ewantRow(1:vecLen) = eosData(eintRad+1:eintRad+vecLen)   ! store desired radiation energy
     call eos_newtonRaphson(vecLen, MODE_DENS_EI_RAD, .TRUE., maskPtr,componentMask)


     eosData(tempRad+1:tempRad+vecLen)=tempRow(1:vecLen)
     if(maskPtr(EOS_PRESRAD)) eosData(presRad+1:presRad+vecLen)=ptotRow(1:vecLen)
     if (eos_forceConstantInput)  then
!!$        eosData(eintRad+1:eintRad+vecLen) = esaveRow(1:vecLen)
     else
        if(maskPtr(EOS_EINTRAD)) eosData(eintRad+1:eintRad+vecLen) = etotRow(1:vecLen)
     end if

     componentMask(:) = (/1,1,0/)
     tempRow(1:vecLen) = tempEleRow(1:vecLen)
     ewantRow(1:vecLen) = eosData(eintEle+1:eintEle+vecLen)+eosData(eintIon+1:eintIon+vecLen)
     call eos_newtonRaphson(vecLen, MODE_DENS_EI_ELE, .TRUE., maskPtr,componentMask)
     eosData(tempEle+1:tempEle+vecLen)=tempRow(1:vecLen)
     if(maskPtr(EOS_PRESELE)) eosData(presEle+1:presEle+vecLen)=ptotRow(1:vecLen)
     if (eos_forceConstantInput)  then
!!$        eosData(eintEle+1:eintEle+vecLen) = esaveRow(1:vecLen)
     else
        if(maskPtr(EOS_EINTELE)) eosData(eintEle+1:eintEle+vecLen) = etotRow(1:vecLen)
     end if

!!$     if(maskPtr(EOS_PRESION)) eosData(presIon+1:presIon+vecLen)=ptotRow(1:vecLen)
     if (eos_forceConstantInput)  then
!!$        eosData(eintIon+1:eintIon+vecLen) = esaveRow(1:vecLen)
     else
        if(maskPtr(EOS_EINTION)) eosData(eintIon+1:eintIon+vecLen) = 0.0
     end if

     if (mode==MODE_DENS_EI_ALL) then
        componentMask(:) = (/1,1,1/)
        tempRow(1:vecLen) = eosData(temp+1:temp+vecLen)
        ewantRow(1:vecLen) = eosData(eint+1:eint+vecLen)   ! re-store desired internal energy
        call eos_newtonRaphson(vecLen, MODE_DENS_EI, .TRUE., maskPtr,componentMask)
     end if
     if (mode==MODE_DENS_EI_GATHER .OR. mode==MODE_DENS_EI_RECAL_GATHER) then
        tempRow(1:vecLen) = eosData(temp+1:temp+vecLen)
        if (ANY(maskPtr((/EOS_PRESELE,EOS_PRESION,EOS_PRESRAD/)) .EQV. .FALSE.)) then
           call Driver_abortFlash("[Eos] cannot calculate MODE_DENS_EI_GATHER without component pressure masks.&
                & Set mask appropriately.")
        end if
        if(.TRUE. .OR. eos_combinedTempRule==0) then ! Force it to get gamc that includes all components.
           ewantRow(1:vecLen) = eosData(eintIon+1:eintIon+vecLen)+ &
                             eosData(eintEle+1:eintEle+vecLen)+ &
                             eosData(eintRad+1:eintRad+vecLen)
           esaveRow(1:vecLen) = ewantRow(1:vecLen)
           call eos_newtonRaphson(vecLen, MODE_DENS_EI, .TRUE., maskPtr,cMask=(/1,1,1/))
#ifdef DEBUG_EOS
           if (any(gamcRow(1:vecLen).LE.1.3)) then
              print*,'gamc:',eosData(gamc+1:gamc+vecLen)
           end if
#endif
           ! DEV: Review this! KW 2016-07-05
           eosData(pres+ilo:pres+ihi) = ( &
             eosData(presEle+ilo:presEle+ihi)+ &
             eosData(presRad+ilo:presRad+ihi))
        else
           eosData(eint+ilo:eint+ihi) = ( &
             eosData(eintEle+ilo:eintEle+ihi)+ &
             eosData(eintRad+ilo:eintRad+ihi))
           
           eosData(pres+ilo:pres+ihi) = ( &
             eosData(presEle+ilo:presEle+ihi)+ &
             eosData(presRad+ilo:presRad+ihi))
        end if
        call setCombinedTemp()
     end if

!!!! NEW EXPERIMENTAL MODE_DENS_EI_MAT_GATHER
!!!!

  case (MODE_DENS_EI_MAT_GATHER)

     eosData(presIon+ilo:presIon+ihi) = 0.0

     eosData(tempIon+ilo:tempIon+ihi) = 0.0

     tempEleRow(1:vecLen) = eosData(tempEle+1:tempEle+vecLen)

     eosData(tempRad+ilo:tempRad+ihi) = sqrt(sqrt( eosData(presRad+ilo:presRad+ihi) / eos_ao3 ))

     componentMask(:) = (/1,1,0/)
     tempRadRow(1:vecLen) = tempEleRow(1:vecLen)
     tempRow(1:vecLen) = tempEleRow(1:vecLen)
     ewantRow(1:vecLen) = eosData(eintEle+1:eintEle+vecLen)+eosData(eintIon+1:eintIon+vecLen)
     call eos_newtonRaphson(vecLen, MODE_DENS_EI_ELE, .TRUE., maskPtr,componentMask)
     eosData(tempEle+1:tempEle+vecLen)=tempRow(1:vecLen)
     if(maskPtr(EOS_PRESELE)) eosData(presEle+1:presEle+vecLen)=ptotRow(1:vecLen)
     if (eos_forceConstantInput)  then
!!$        eosData(eintEle+1:eintEle+vecLen) = esaveRow(1:vecLen)
     else
        if(maskPtr(EOS_EINTELE)) eosData(eintEle+1:eintEle+vecLen) = etotRow(1:vecLen)
     end if

!!$     if(maskPtr(EOS_PRESION)) eosData(presIon+1:presIon+vecLen)=ptotRow(1:vecLen)
     if (eos_forceConstantInput)  then
!!$        eosData(eintIon+1:eintIon+vecLen) = esaveRow(1:vecLen)
     else
        if(maskPtr(EOS_EINTION)) eosData(eintIon+1:eintIon+vecLen) = 0.0
     end if

!!     eosData(entr+ilo:entr+ihi) = 0.0

     if (.TRUE.) then
        if (ANY(maskPtr((/EOS_PRESELE,EOS_PRESION,EOS_PRESRAD/)) .EQV. .FALSE.)) then
           call Driver_abortFlash("[eos_mgamma] cannot calculate MODE_DENS_EI_GATHER without component pressure masks.&
                & Set mask appropriately.")
        end if

        eosData(eint+ilo:eint+ihi) = ( &
             eosData(eintEle+ilo:eintEle+ihi)+ &
             eosData(eintRad+ilo:eintRad+ihi))
        eosData(presRad+ilo:presRad+ihi) = eosData(eintRad+ilo:eintRad+ihi) * eosData(dens+ilo:dens+ihi) / 3.0
        eosData(tempRad+ilo:tempRad+ihi) = sqrt(sqrt( eosData(presRad+ilo:presRad+ihi) / eos_ao3 ))

        eosData(pres+ilo:pres+ihi) = ( &
             eosData(presEle+ilo:presEle+ihi)+ &
             eosData(presRad+ilo:presRad+ihi))
        call setCombinedTemp()
     end if

    !! Note that we make here a choice for the arbitrary additive constant in the electron entropy 
!!$     if(maskPtr(EOS_ENTRELE)) &
!!$          call setEleEntropy(entrele,tempEle)
!!$     if(maskPtr(EOS_ENTRRAD)) then
!!$        call setRadEntropy(entrrad,tempRad,eintRad)
!!$     end if

!!!! *ADDITIONAL* NEW EXPERIMENTAL MODE_DENS_EI_MAT_EQUI
!!!!

  case (MODE_DENS_EI_MAT_EQUI)

     eosData(presIon+ilo:presIon+ihi) = 0.0

     eosData(tempIon+ilo:tempIon+ihi) = 0.0

     eosData(tempRad+ilo:tempRad+ihi) = sqrt(sqrt( eosData(presRad+ilo:presRad+ihi) / eos_ao3 ))


!!     eosData(entr+ilo:entr+ihi) = 0.0



     if (.TRUE.) then
        if (ANY(maskPtr((/EOS_PRESELE,EOS_PRESION,EOS_PRESRAD/)) .EQV. .FALSE.)) then
           call Driver_abortFlash("[eos_mgamma] cannot calculate MODE_DENS_EI_GATHER without component pressure masks.&
                & Set mask appropriately.")
        end if
        eosData(presRad+ilo:presRad+ihi) = 0.0
        eosData(eint+ilo:eint+ihi) = ( &
             eosData(eintEle+ilo:eintEle+ihi))

        call eos_newtonRaphson(vecLen, MODE_DENS_EI, .TRUE., maskPtr,cMask=(/1,1,0/))
        eosData(entr+ilo:entr+ihi)=stotRow(1:rowLen) !(replace with expression based on component states)
        eosData(eint+ilo:eint+ihi) = ( &
             eosData(eintEle+ilo:eintEle+ihi)+ &
             eosData(eintRad+ilo:eintRad+ihi))
        eosData(presRad+ilo:presRad+ihi) = eosData(eintRad+ilo:eintRad+ihi) * eosData(dens+ilo:dens+ihi) / 3.0
        eosData(tempRad+ilo:tempRad+ihi) = sqrt(sqrt( eosData(presRad+ilo:presRad+ihi) / eos_ao3 ))

        if(eos_combinedTempRule==0) then
        endif

        eosData(pres+ilo:pres+ihi) = ( &
             eosData(presEle+ilo:presEle+ihi)+ &
             eosData(presRad+ilo:presRad+ihi))
        call setCombinedTemp()
     end if

    !! Note that we make here a choice for the arbitrary additive constant in the electron entropy 
!!$     if(maskPtr(EOS_ENTRELE)) &
!!$          call setEleEntropy(entrele,tempEle)
!!$     if(maskPtr(EOS_ENTRRAD)) then
!!$        call setRadEntropy(entrrad,tempRad,eintRad)
!!$     end if



!!!!
!!!! END OF NEW EXPERIMENTAL MODE_DENS_EI_MAT

     !==============================================================================

     !      MODE_DENS_PRES  pressure and density given

  case(MODE_DENS_PRES)

     componentMask(:) = 1
     call eos_newtonRaphson(vecLen, mode, .FALSE., maskPtr, componentMask)




     !==============================================================================

     ! Unknown EOS mode selected

  case default
     if (eos_meshMe .EQ. MASTER_PE) print*, '[eos_helmholtz] Error: unknown input mode', mode
     call Driver_abortFlash('[Eos] Error: unknown input mode in subroutine eos_helmholtz')
  end select


  ! Get the optional values
  if(present(mask)) then
     ! Entropy derivatives
     if(mask(EOS_DST)) then
        dst = (EOS_DST-1)*vecLen
        eosData(dst+1:dst+vecLen) = dstRow(1:vecLen)
     end if
     if(mask(EOS_DSD)) then
        dsd = (EOS_DSD-1)*vecLen
        eosData(dsd+1:dsd+vecLen) = dsdRow(1:vecLen)
     end if
     if(mask(EOS_DPT)) then
        dpt = (EOS_DPT-1)*vecLen
        eosData(dpt+1:dpt+vecLen) = dptRow(1:vecLen)
     end if
     if(mask(EOS_DPD)) then
        dpd = (EOS_DPD-1)*vecLen
        eosData(dpd+1:dpd+vecLen) = dpdRow(1:vecLen)
     end if
     if(mask(EOS_DET))then
        det = (EOS_DET-1)*vecLen
        eosData(det+1:det+vecLen) = detRow(1:vecLen)
     end if
     if(mask(EOS_DED))then 
        ded = (EOS_DED-1)*vecLen
        eosData(ded+1:ded+vecLen) = dedRow(1:vecLen)
     end if
     if(mask(EOS_DEA))then 
        dea = (EOS_DEA-1)*vecLen
        eosData(dea+1:dea+vecLen) = deaRow(1:vecLen)
     end if
     if(mask(EOS_DEZ))then 
        dez = (EOS_DEZ-1)*vecLen
        eosData(dez+1:dez+vecLen) = dezRow(1:vecLen)
     end if
     if(mask(EOS_PEL))then 
        pel = (EOS_PEL-1)*vecLen
        eosData(pel+1:pel+vecLen) = pelRow(1:vecLen)
     end if
     if(mask(EOS_NE))then 
        ne = (EOS_NE-1)*vecLen
        eosData(ne+1:ne+vecLen) = neRow(1:vecLen)
     end if
     if(mask(EOS_ETA))then 
        eta = (EOS_ETA-1)*vecLen
        eosData(eta+1:eta+vecLen) = etaRow(1:vecLen)
     end if

     if(mask(EOS_CV))then
        if(mask(EOS_DET)) then
           c_v = (EOS_CV-1)*vecLen
           eosData(c_v+1:c_v+vecLen) = cvRow(1:vecLen)
        else
           call Driver_abortFlash("[eos_helmholtz] cannot calculate C_V without DET.  Set mask appropriately.")
        end if
     end if

     if(mask(EOS_CP))then
        if(mask(EOS_CV).and.mask(EOS_DET)) then
           c_p = (EOS_CP-1)*vecLen
           eosData(c_p+1:c_p+vecLen) = cpRow(1:vecLen)
        else
           call Driver_abortFlash("[eos_helmholtz] cannot calculate C_P without C_V and DET.  Set mask appropriately.")
        end if
     end if
  end if

  !! DEV: Check whether cvIonRow and cvEleRow are computed correctly in eos_helm.F90.

    if(mask(EOS_CVELE))then
       c_v = (EOS_CVELE-1)*vecLen
       eosData(c_v+1:c_v+vecLen) = cvEleRow(1:vecLen)
     end if

    if(mask(EOS_CVION))then
       c_v = (EOS_CVION-1)*vecLen
       eosData(c_v+1:c_v+vecLen) = cvIonRow(1:vecLen)
     end if
  !!


  if (arrayBoundHackMode .eqv. .true.) then
     deallocate(maskPtr)
  end if


  !! Close up arrays if previously allocated
#ifndef FIXEDBLOCKSIZE  
  call eos_vecDealloc()
#endif

  return

!! MANOS B
contains

#if(0)
  subroutine setCombinedGamc(gamcCombined)
    integer, intent(in) :: gamcCombined

    eosData(gamcCombined+ilo:gamcCombined+ihi) = ( &
        gamIon(ilo:ihi)*eosData(presEle+ilo:presEle+ihi)/Zp + &
           eos_gammaEle*eosData(presEle+ilo:presEle+ihi)*dynamicZ/Zp + &
           eos_gammaRad*eosData(presRad+ilo:presRad+ihi) ) / eosData(pres+ilo:pres+ihi)
!!$    print*,'setCGamc->',eosData(gamcCombined+ilo:gamcCombined+ihi)

  end subroutine setCombinedGamc

  subroutine setCombinedGamc1(gamcCombined,ioff)
    integer, intent(in) :: gamcCombined,ioff

    eosData(gamcCombined+ioff) = ( &
        ( gamIon(ioff) + eos_gammaEle*eosdata(zbar+ioff) )&
          * eosData(presEle+ioff)/(eosData(zbar+ioff)+1.0) + &
           eos_gammaRad*eosData(presRad+ioff) ) / eosData(pres+ioff)
!!$    print*,'setCGamc->',eosData(gamcCombined+ilo:gamcCombined+ihi)

  end subroutine setCombinedGamc1
#endif

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

  subroutine eos_newtonRaphson(vecLen, mode, enerWanted, eosMask,cMask)
    implicit none
    integer,intent(IN) :: vecLen, mode
    logical,intent(IN) :: enerWanted
    logical,optional, dimension(EOS_VARS+1:EOS_NUM),INTENT(in)::eosMask
    integer,optional, dimension(EOSCOMP_NUM_COMPONENTS),INTENT(in)::cMask

    real, pointer, dimension(:) :: xXtotRow, xXwantRow
    real    :: told
    real,dimension(vecLen):: bnd_lo, bnd_hi

    if (enerWanted) then
       !  ewantRow is our desired EI input
       xXtotRow  => etotRow
       xXwantRow => ewantRow
    else                    !pressure is wanted
       xXtotRow  => ptotRow
       xXwantRow => pwantRow
    end if

    ! Initialize the errors
    error(:) = 0.0e0

    !! set wide bracket bounds
    bnd_lo(:) = eos_smallt
    bnd_hi(:) = eos_largeT

    ! Do the first eos call with all the zones in the pipe
    !  NOTE that eos_helm can ONLY operate in the equivalent of
    !  MODE_DENS_TEMP, as it returns pressure, energy and derivatives only
    !  So if you send in a crappy temperature here, you'll get a crappy starting
    !  position and the iteration won't converge.
    !  Initial temperature here is what is stored in the grid, even though we 
    !    SUSPECT this is not in equilibrium (or we wouldn't be calling Eos if it was fine)
    if (.NOT. all(tempRow(1:vecLen) > 0)) &
         print*,'eos_nr WARN:',tempRow(1:vecLen)
    ! Catch really terrible guesses.
    do k = vecBegin, vecEnd
       if (tempRow(k) > eos_largeT) tempRow(k) = eos_largeT
       if (tempRow(k) < eos_smallt) tempRow(k) = eos_smallt
    enddo
    if (.NOT. present(cMask)) then
       tempRadRow(1:vecLen) = tempRow(1:vecLen)
       tempIonRow(1:vecLen) = tempRow(1:vecLen)
       tempEleRow(1:vecLen) = tempRow(1:vecLen)
    else
       if (cMask(EOSCOMP_RAD).NE.0) tempRadRow(1:vecLen) = tempRow(1:vecLen)
       if (cMask(EOSCOMP_ION).NE.0) tempIonRow(1:vecLen) = tempRow(1:vecLen)
       if (cMask(EOSCOMP_ELE).NE.0) tempEleRow(1:vecLen) = tempRow(1:vecLen)
    end if

    call eos_helm(vecBegin,vecEnd,eosMask,componentMask=cMask)
    !  Now eos_helm has returned ptotRow, etotRow, dXXtRow, and gamcRow


    !  Create initial condition
    do k = vecBegin, vecEnd
       if (enerWanted) then
          !  ewantRow is our desired EI input
          tnew(k) = tempRow(k) - (etotRow(k) - ewantRow(k)) / detRow(k)
       else                    !pressure is wanted
          tnew(k) = tempRow(k) - (ptotRow(k) - pwantRow(k)) / dptRow(k)
       end if

       ! if we jump out of the brackets, check if there is a solution
       if ( tnew(k) < bnd_lo(k) ) then
          told = tempRow(k)
          tempRow(k) = bnd_lo(k)
          call eos_helm(k,k,eosMask,componentMask=cMask)
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
          call eos_helm(k,k,eosMask,componentMask=cMask)
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
       ! (This MAY not be needed anymore with bracketing, byt MAY still reduce
       !  the number of iterations in some bad cases (?) )
       if (tnew(k) .GT. 10.e0*tempRow(k)) tnew(k) =  & 
            &           10.e0*tempRow(k)
       if (tnew(k) .LT. 0.1e0*tempRow(k)) tnew(k) =  & 
            &           0.1e0*tempRow(k)

       ! Compute the error
       error(k) = abs((tnew(k) - tempRow(k)) / tempRow(k))

       ! Store the new temperature
       ! but not if the guess was already good enough
       ! this prevents values from shifting if the EOS is called more than once
       ! without any changes in between (e.g. in guard cells)
       if ( .not. (error(k)< eos_tol) ) then
          tempRow(k) = tnew(k)
       endif

       ! DEV: Disable the following 'freezing' logic? - KW
       ! Check if we are freezing, if so set the temperature to smallt, and adjust 
       ! the error so we don't wait for this one
       if (tempRow(k) .LT. eos_smallt) then
          tempRow(k) = eos_smallt
          error(k)    = 0.1*eos_tol
       endif

    enddo

    ! Loop over the zones individually now
    do k = vecBegin, vecEnd
       do i = 2, eos_maxNewton
          if (error(k) < eos_tol) goto 70

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
          call eos_helm(k,k,eosMask,componentMask=cMask)

          if (enerWanted) then
             tnew(k) = tempRow(k) - (etotRow(k) - ewantRow(k)) / detRow(k)
          else                    !pressure is wanted
             tnew(k) = tempRow(k) - (ptotRow(k) - pwantRow(k)) / dptRow(k)
          end if
          ! if we jump out of the brackets, check if there is a solution
          if ( tnew(k) < bnd_lo(k) ) then
             told = tempRow(k)
             tempRow(k) = bnd_lo(k)
             call eos_helm(k,k,eosMask,componentMask=cMask)
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
          else if ( tnew(k) > bnd_hi(k) ) then
             told = tempRow(k)
             tempRow(k) = bnd_hi(k)
             call eos_helm(k,k,eosMask,componentMask=cMask)
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
          ! (this is probably not needed anymore with bracketing, but might reduce
          !  the number of iterations in some bad cases (?) )
          if (tnew(k) .GT. 10.e0*tempRow(k)) tnew(k) =  & 
               &              10.e0*tempRow(k)
          if (tnew(k) .LT. 0.1e0*tempRow(k)) tnew(k) =  & 
               &              0.1e0*tempRow(k)

          ! Compute the error
          error(k) = abs((tnew(k) - tempRow(k)) / tempRow(k))

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
       print *, 'Newton-Raphson failed in subroutine eos_helmholtz'

       select case(mode)
       case(MODE_DENS_EI)
          print *, '(e and rho as input):'
          print *, ' ewant= ', ewantRow(k)
       case(MODE_DENS_PRES)
          print *, '(pres and rho as input):'
          print *, ' pwant= ', pwantRow(k)

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

       case default
          print *, 'Unusual Eos mode:', mode
       end select
       print *, ' '
       print *, 'too many iterations', eos_maxNewton
       print *, ' '
       print *, ' k    = ', k,vecBegin,vecEnd
       print *, ' temp = ', tempRow(k)
       print *, ' dens = ', denRow(k)
       print *, ' abar = ', abarRow(k)
       print *, ' zbar = ', zbarRow(k)
       print *, ' pres = ', ptotRow(k)

       call Driver_abortFlash('[Eos] Error: too many iterations in Newton-Raphson')


       ! Land here if the Newton iteration converged
       !  jumps out of the iterations, but then continues to the next vector location

70       continue           
!   debugging checks for reporting hitting boundaries disabled by default
#ifdef DEBUG_EOS
       if ( tempRow(k) == eos_smallt) then
          print *, 'eos: at dens =', denRow(k), ' shooting for e = ', ewantRow(k), &
                  'needed faking temperature at smallt= ', tempRow(k)
       endif
       if (bnd_hi(k) == eos_largeT) then
          if ( abs(tempRow(k) - eos_largeT)/eos_largeT < 2.0*eos_tol ) then
             print *, 'eos: at dens =', denRow(k), ' shooting for e = ', ewantRow(k), &
                  'needed faking temperature at largeT= ', tempRow(k)
          endif
       endif
#endif

    end do

    ! Crank through the entire eos one last time

    if (.NOT. present(cMask)) then
       tempRadRow(1:vecLen) = tempRow(1:vecLen)
       tempIonRow(1:vecLen) = tempRow(1:vecLen)
       tempEleRow(1:vecLen) = tempRow(1:vecLen)
    else
       if (cMask(EOSCOMP_RAD).NE.0) tempRadRow(1:vecLen) = tempRow(1:vecLen)
       if (cMask(EOSCOMP_ION).NE.0) tempIonRow(1:vecLen) = tempRow(1:vecLen)
       if (cMask(EOSCOMP_ELE).NE.0) tempEleRow(1:vecLen) = tempRow(1:vecLen)
    end if
    call eos_helm(vecBegin,vecEnd,eosMask,componentMask=cMask)

    ! Fill the FLASH arrays with the results.  
    !  In MODE_DENS_EI, we should be generating temperature and pressure (plus gamma and entropy)

    select case(mode)
    case(MODE_DENS_EI)
       eosData(temp+1:temp+vecLen)=tempRow(1:vecLen)
       eosData(pres+1:pres+vecLen)=ptotRow(1:vecLen)
       !  Update the energy to be the true energy, instead of the energy we were trying to meet
       !  ConstantInput LBR and KW believe this is WRONG -- the input arrays should not be changed
       if (eos_forceConstantInput)  then
          eosData(eint+1:eint+vecLen) = esaveRow(1:vecLen)
       else
          eosData(eint+1:eint+vecLen) = etotRow(1:vecLen)
       end if
    case(MODE_DENS_PRES)
       eosData(temp+1:temp+vecLen)=tempRow(1:vecLen)
       ! Update the pressure to be the equilibrium pressure, instead of the pressure we were trying to meet
       !  ConstantInput LBR and KW believe this is wrong.  See notes at the top of the routine
       if (eos_forceConstantInput) then
          eosData(pres+1:pres+vecLen) = psaveRow(1:vecLen)
       else
          eosData(pres+1:pres+vecLen) = ptotRow(1:vecLen)
       end if
       eosData(eint+1:eint+vecLen)=etotRow(1:vecLen)
    end select

    eosData(gamc+1:gamc+vecLen)=gamcRow(1:vecLen)
    eosData(entr+1:entr+vecLen)=stotRow(1:vecLen)

  end subroutine eos_newtonRaphson
!! MANOS E

end subroutine eos_helmholtz


