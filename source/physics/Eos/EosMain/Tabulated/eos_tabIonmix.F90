!!****if* source/physics/Eos/EosMain/Tabulated/eos_tabIonmix
!!
!! NAME
!!
!!  eos_tabIonmix
!!
!! SYNOPSIS
!!
!!  call eos_tabIonmix( integer(IN) :: mode,
!!                      integer(IN) :: vecLen,
!!                      real(INOUT) :: eosData(EOS_NUM*vecLen),
!!            optional, integer(IN) :: vecBegin,
!!            optional, integer(IN) :: vecEnd,
!!            optional, integer(IN) :: eosType,
!!            optional, integer(IN) :: subtype,
!!            optional, integer(IN) :: material,
!!      optional,target,logical(IN) :: mask(EOS_VARS+1:EOS_NUM))
!!                 
!!
!! DESCRIPTION
!!
!!  This routine implements the tabulated version of the equation of state
!!  for IONMIX tables.
!!
!!  ARGUMENTS
!!
!!  mode :    Selects the mode of operation of the Eos unit.
!!             The valid values are MODE_DENS_EI, MODE_DENS_PRES and  
!!             MODE_DENS_TEMP as decribed above.
!!
!!  vecLen   : number of points for each input variable
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
!!
!!  vecEnd   : Index of last cell in eosData to handle. 
!!             Can be used to limit operation to a subrange of cells, untested.
!!             If not present, the default is vecLen.
!!
!!  eosType  : the type of eos to be applied
!!
!!  material : Indicates to which material the EOS is to be applied, 
!!             in a multi-type multi-material context.
!!             Given as an index into the UNK solution vector,
!!             if valid we should have SPECIES_BEGIN <= material <= SPECIES_END.
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
!!
!!
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

subroutine eos_tabIonmix(mode, vecLen, eosData, vecBegin, vecEnd, eosType, subtype, material, mask)

!==============================================================================
  use Eos_data, ONLY : eos_gasConstant, &
       eos_meshMe
  use eos_helmConstData, ONLY: eos_ao3, eos_avo
  use eos_tabInterface, ONLY: eos_tabGetSpeciesTabulatedData, eos_tabIonmix4
  use eos_tabData, ONLY: EOS_TAB_NCOMP,EOS_TAB_NALLTAB, &
                          EOS_TAB_FOR_ION, EOS_TAB_FOR_ELE, EOS_TAB_FOR_MAT, &
                          EOS_TABVT_ZF, EOS_TABVT_EN, EOS_TABVT_PR, &
                          EOS_TABVT_ENTR, &
                          EOS_TABINT_DERIV_0, EOS_TAB_NDERIVS, &
                          EOS_TABINT_DERIV_DT, EOS_TABINT_DERIV_DD

  implicit none
#include "constants.h"
#include "Eos.h"
#include "Flash.h"
#include "Multispecies.h"

  !     Arguments
  integer, INTENT(in) :: mode, vecLen
  real,INTENT(inout), dimension(EOS_NUM*vecLen) :: eosData 
  integer,optional,INTENT(in) :: vecBegin,vecEnd
  integer,optional,INTENT(in) :: eosType,subtype
  integer,optional,INTENT(in) :: material
  logical, optional, INTENT(in),target,dimension(EOS_VARS+1:EOS_NUM) :: mask

  ! This is the variable that is used internally -- set to false unless mask comes in
  logical,save,target, dimension(EOS_VARS+1:EOS_NUM) :: maskInternal
  data maskInternal / EOS_DERIVS * .FALSE. /
  logical,pointer, dimension(:) :: maskPtr

  real,dimension(vecLen) :: gamIon, gamM1Ion, ggprodIon
  integer :: dens, temp, pres, eint, abar, zbar
  integer :: entr, dst, dsd
  integer :: dpt, dpd, det, ded, c_v, c_p, gamc, game, pel, ne, eta
  integer :: cvion, cvele
  integer :: tempIon,tempEle,tempRad                                                
  integer :: eintIon,eintEle,eintRad                                                
  integer :: presIon,presEle,presRad
  integer :: entrEle
  integer :: i, ilo,ihi, rowLen

  logical :: arrayBoundHackMode
  integer :: myType, mySubtype, species, specno
  logical :: useNRForAddtlCombined !use vectors set by NR call for addtl. combined components
  !addtl. comp. = cp, cv, ...
  logical :: useNRForGamcCombined  !use vectors set by NR call for gamc combined components
  logical :: doUncoupledForAddtl ! assume components thermally uncoupled for addtl. out variables
  logical :: doUncoupledForGamc  ! assume components thermally uncoupled for gamc out variables
  logical :: doCoupledForGamc    ! handle components as thermally coupled for gamc out variables
  logical :: wantComb,wantIon,wantEle,wantRad
  integer :: tempToUse, combTable

  integer :: maxDerivsE, maxDerivsZ
  logical :: wanted(EOS_TAB_NCOMP,EOS_TAB_NALLTAB)
  integer :: neededTabDerivs(EOS_TAB_NCOMP,EOS_TAB_NALLTAB)
  real :: tabData(0:EOS_TAB_NDERIVS,1:EOS_TABVT_ENTR,EOS_TAB_NCOMP)

  maskPtr => maskInternal
  if (present(mask)) then
     maskPtr => mask
  end if
 

  if (present(eosType)) then
     myType = eosType
  else
     myType = EOS_TAB
  end if
  select case (myType)
     case(EOS_TAB)
        ! ok...
     case default
        if (eos_meshMe==MASTER_PE) print*,"eos_tabIonmix: only supports EOS type EOS_TAB, called for eosType=", myType
        call Driver_abortFlash("eos_tabIonmix: only supports EOS_TAB, called for wrong EOS type!")
  end select

  if (present(subtype)) then
     mySubtype = subtype
  else
     mySubtype = 1
  end if
  select case (mySubtype)
     case(1)
        ! ok, continue below...
     case(4,6)
        call eos_tabIonmix4(mode, vecLen, eosData, vecBegin, vecEnd, eosType, subtype, material, mask)
        return                  ! work done, RETURN IMMEDIATELY
     case default
        if (eos_meshMe==MASTER_PE) print*,"eos_tabIonmix: only supports EOS subtypes 1 and 4, called for subtype=", mySubtype
        call Driver_abortFlash("eos_tabIonmix: only supports EOS subtypes 1 and 4, called for wrong EOS subtype!")
  end select

  if (lbound(maskPtr,1) == 1) then
     !Some older compiler versions do not set the lower bound of maskPtr to EOS_VARS+1.
     !This is wrong, but we work around it here.
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

  species = 0
  if (present(material)) then
     if (NSPECIES < 1) then
        species = -2
     else if (material < SPECIES_BEGIN) then
        species = -3
     else if (material > SPECIES_END) then
        species = -1
     else
        species = material
     end if
  end if
  if (species .LE. 0) then
     if (eos_meshMe==MASTER_PE) print*,"eos_tabIonmix: invalid species=",species
     call Driver_abortFlash("eos_tabIonmix: species number not given or invalid!")
  end if

  specno = species - SPECIES_BEGIN + 1


#ifdef DEBUG_EOS
  print*,'Calling eos_tabIonmix, mode',mode,' specno',specno
#endif

  useNRForAddtlCombined = .FALSE.; useNRForGamcCombined = .FALSE.
  doUncoupledForAddtl = .FALSE.; doUncoupledForGamc = .FALSE.; doCoupledForGamc = .FALSE.

#define dynamicZ eosData(zbar+ilo:zbar+ihi)
#define Zp (dynamicZ+1)

!============================================================================

  pres = (EOS_PRES-1)*vecLen
  dens = (EOS_DENS-1)*vecLen
  temp = (EOS_TEMP-1)*vecLen
  eint = (EOS_EINT-1)*vecLen
  gamc = (EOS_GAMC-1)*vecLen
  game = gamc                   !for now
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
  det = (EOS_DET-1)*vecLen
  ded = (EOS_DED-1)*vecLen
  dpt = (EOS_DPT-1)*vecLen
  dpd = (EOS_DPD-1)*vecLen
  cvion = (EOS_CVION-1)*vecLen
  cvele = (EOS_CVELE-1)*vecLen


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
     print*,'[eos_tabIonmix3T] ilo is',ilo
     call Driver_abortFlash("[eos_tabIonmix3T] invalid ilo")
  end if
  if (ihi < 1 .OR. ihi > vecLen) then
     print*,'[eos_tabIonmix3T] ihi is',ihi
     call Driver_abortFlash("[eos_tabIonmix3T] invalid ihi")
  end if
  if (rowLen < 0 .OR. rowLen > vecLen) then
     print*,'[eos_tabIonmix3T] rowLen is',rowLen
     call Driver_abortFlash("[eos_tabIonmix3T] invalid rowLen")
  end if
#endif
  if (rowLen == 0) then
     print*,'[eos_tabIonmix3T] rowLen is 0.'
  end if

  select case (mode)
  case(MODE_DENS_TEMP_ALL)
  case(MODE_DENS_TEMP_COMP,MODE_DENS_TEMP_GATHER)
     useNRForAddtlCombined = .FALSE.; useNRForGamcCombined = .FALSE.
  case(MODE_DENS_TEMP,MODE_DENS_TEMP_EQUI)
!!$     useNRForAddtlCombined = .TRUE.; useNRForGamcCombined = .TRUE. !NONO in this file!
     doCoupledForGamc = .TRUE.

  case(MODE_DENS_EI_ALL)
  case(MODE_DENS_EI_COMP, MODE_DENS_EI_GATHER, &
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

  
  gamIon(ilo:ihi) = 5./3.
  gamM1Ion(ilo:ihi) = 1.5
  ggprodIon(ilo:ihi) = gamM1Ion(ilo:ihi) * eos_gasConstant

#if 0
  if(maskPtr(EOS_ENTRELE)) then
     call Driver_abortFlash('eos_tabIonmix: computing of electron entropy is not supported!')
  end if
#endif

!!$  wanted(:,:) = .FALSE.
!!$  wanted(EOS_TAB_FOR_ELE,EOS_TABVT_ZF) = .TRUE.
  neededTabDerivs(:,:) = -1        !By default we do not even need "zeroth deriv", i.w., table function itself.
  neededTabDerivs(EOS_TAB_FOR_ELE,EOS_TABVT_ZF) = 0

  !Is MODE_DENS_TEMP_GATHER to be considered invalid here?
  !Is MODE_DENS_TEMP_COMP to be considered invalid here? (because more than one temp involved?)
  !Is MODE_DENS_TEMP_ALL to be considered invalid here? (because more than one temp involved?)
  if (mode==MODE_DENS_TEMP .OR. mode==MODE_DENS_TEMP_ALL .OR. &
       mode==MODE_DENS_TEMP_MAT_EQUI .OR. &
       mode==MODE_DENS_TEMP_EQUI) then
     wantComb = .TRUE.
!!$     wanted(EOS_TAB_FOR_MAT,EOS_TABVT_EN) = .TRUE.
     neededTabDerivs(EOS_TAB_FOR_MAT,EOS_TABVT_EN) = 0
  else
     wantComb = .FALSE.
  end if
!!$  if (mode==MODE_DENS_TEMP_ION .OR. mode==MODE_DENS_TEMP_ALL .OR. &
!!$       mode==MODE_DENS_TEMP_MAT_EQUI .OR. &
!!$       mode==MODE_DENS_TEMP_COMP .OR. mode==MODE_DENS_TEMP_EQUI) then
!!$     wantIon = .TRUE.
!!$  else
     wantIon = .FALSE.          !No, do not want ion data *from table*
!!$  end if
  if (mode==MODE_DENS_TEMP_ELE .OR. mode==MODE_DENS_TEMP_ALL .OR. &
       mode==MODE_DENS_TEMP_MAT_EQUI .OR. &
       mode==MODE_DENS_TEMP_COMP .OR. mode==MODE_DENS_TEMP_EQUI) then
     wantEle = .TRUE.
!!$     wanted(EOS_TAB_FOR_MAT,EOS_TABVT_EN) = .TRUE.
     neededTabDerivs(EOS_TAB_FOR_MAT,EOS_TABVT_EN) = 0
  else
     wantEle = .FALSE.
  end if
  if (mode==MODE_DENS_TEMP_RAD .OR. mode==MODE_DENS_TEMP_ALL .OR. &
       mode==MODE_DENS_TEMP_COMP .OR. mode==MODE_DENS_TEMP_EQUI) then
     wantRad = .TRUE.
  else
     wantRad = .FALSE.
  end if
  if (mode==MODE_DENS_TEMP_ION) then
     tempToUse = tempIon
     combTable = EOS_TAB_FOR_ION !DEV: combTable probably obsolete in this file
  else if (mode==MODE_DENS_TEMP_ELE) then
     tempToUse = tempEle
     combTable = EOS_TAB_FOR_ELE
  else if (mode==MODE_DENS_TEMP_RAD) then
     tempToUse = tempRad
     combTable = EOS_TAB_FOR_ELE !???
  else
     tempToUse = temp
     combTable = EOS_TAB_FOR_MAT
  end if

  if(mode==MODE_DENS_TEMP_EQUI) then
     !Equalize temperature.  Currently done below.
  end if

  maxDerivsE = 0
  if(maskPtr(EOS_DET).OR. .TRUE.) then
     maxDerivsE = 1
     neededTabDerivs(EOS_TAB_FOR_MAT,EOS_TABVT_EN) = 1
  end if
  if(maskPtr(EOS_DED)) then
     maxDerivsE = 2
     neededTabDerivs(EOS_TAB_FOR_MAT,EOS_TABVT_EN) = 2
  end if
  maxDerivsZ = 0
  if(wantComb .OR. wantEle) then
     maxDerivsZ = 2
     neededTabDerivs(EOS_TAB_FOR_ELE,EOS_TABVT_ZF) = 2
  end if

  do i=ilo,ihi
     call eos_tabGetSpeciesTabulatedData(species=specno, &
          speciesTemperature=eosData(tempToUse+i),&
          speciesDensity=eosData(dens+i)/eosData(abar+i)*eos_avo, &
          wantedDerivs=neededTabDerivs,&
!!$          needZFDerivs=maxDerivsZ,&
!!$          needENDerivs=maxDerivsE,&
!!$          needHCDerivs=-1,&
          outData=tabData)
     where (transpose(neededTabDerivs(:,1:EOS_TABVT_ENTR)).GE.2)
        tabData(2,:,:) = tabData(2,:,:)/eosData(abar+i)*eos_avo
     end where

     eosData(zbar+i) = max(eosData(zbar+i), &
                           tabData(EOS_TABINT_DERIV_0,EOS_TABVT_ZF,EOS_TAB_FOR_ELE))
     if (eosData(zbar+i) == 0.0) print*,'zbar is 0.',EOS_TAB_FOR_ION,EOS_TAB_FOR_ELE,EOS_TAB_FOR_MAT
     if (wantComb) then
        eosData(eint+i) = tabData(EOS_TABINT_DERIV_0,EOS_TABVT_EN,EOS_TAB_FOR_MAT) + &
             eosData(eintRad+i)
     end if
     if (wantIon) then
        if(maskPtr(EOS_EINTION)) &
             eosData(eintIon+i) = tabData(EOS_TABINT_DERIV_0,EOS_TABVT_EN,EOS_TAB_FOR_ION)
     end if
     if (wantEle) then
        if(maskPtr(EOS_EINTELE)) then
           eosData(eintEle+i) = tabData(EOS_TABINT_DERIV_0,EOS_TABVT_EN,EOS_TAB_FOR_MAT)
           eosData(eintEle+i) = eosData(eintEle+i) - ggprodIon(i) * eosData(tempToUse+i) / eosData(abar+i)
        end if
     end if
#ifdef RAD_HERE
     if (wantRad) then
        if(maskPtr(EOS_EINTRAD)) then
           eosData(eintRad+i) = &
                3.0e0 * eos_ao3 * eosData(tempToUse+i)**4 / eosData(dens+i)
        end if
     end if
#endif
     if(maskPtr(EOS_CVION)) then
        if (wantComb) then
           eosData(cvion+i) = ggprodIon(i) / eosData(abar+i)
        else if (wantEle) then
           eosData(cvion+i) = 0.0
        else if (wantRad) then
           eosData(cvion+i) = 0.0
        else                    !wants ions, I guess!
           eosData(cvion+i) = ggprodIon(i) / eosData(abar+i)
        end if
     end if
     if(maskPtr(EOS_CVELE)) then
        if (wantComb .OR. wantEle) then
           eosData(cvele+i) = tabData(EOS_TABINT_DERIV_DT,EOS_TABVT_EN,EOS_TAB_FOR_MAT)
           eosData(cvele+i) = eosData(cvele+i) - ggprodIon(i) / eosData(abar+i)
        else if (wantRad) then
           eosData(cvele+i) = 0.0
        else                    !wants ions, I guess!
           eosData(cvele+i) = 0.0
        end if
     end if
     if(maskPtr(EOS_DET).or. .TRUE.) then
        if (wantComb) then
           eosData(det+i) = tabData(EOS_TABINT_DERIV_DT,EOS_TABVT_EN,EOS_TAB_FOR_MAT)
#ifdef RAD_HERE
           eosData(det+i) = &
                eosData(det+i) +  12.0e0 * eos_ao3 * eosData(tempToUse+i)**3 / eosData(dens+i)
#else
           eosData(det+i) = &
                eosData(det+i) +  4 * eosData(eintRad+i) / eosData(tempToUse+i)
#endif
        else if (wantEle) then
           eosData(det+i) = tabData(EOS_TABINT_DERIV_DT,EOS_TABVT_EN,EOS_TAB_FOR_MAT)
           eosData(det+i) = eosData(det+i) - ggprodIon(i) / eosData(abar+i)
        else if (wantRad) then
#ifdef RAD_HERE
           eosData(det+i) = 12.0e0 * eos_ao3 * eosData(tempToUse+i)**3 / eosData(dens+i)
#else
           eosData(det+i) = 4 * eosData(eintRad+i) / eosData(tempToUse+i)
#endif
        else                    !wants ions, I guess!
           eosData(det+i) = ggprodIon(i) / eosData(abar+i)
        end if
     end if
     if(maskPtr(EOS_DED)) then
        if (wantComb) then
           eosData(ded+i) = tabData(EOS_TABINT_DERIV_DD,EOS_TABVT_EN,EOS_TAB_FOR_MAT)
#ifdef RAD_HERE
           eosData(ded+i) = eosData(ded+i) - 3.0e0 * eos_ao3 * &
                (eosData(tempToUse+i)*eosData(tempToUse+i) / eosData(dens+i))**2
#else
           eosData(ded+i) = eosData(ded+i) - eosData(eintRad+i) / eosData(dens+i) ! DEV: THIS dens?? - KW
#endif
        else if (wantEle) then
           eosData(ded+i) = tabData(EOS_TABINT_DERIV_DD,EOS_TABVT_EN,EOS_TAB_FOR_MAT)
        else if (wantRad) then
#ifdef RAD_HERE
           eosData(ded+i) = - 3.0e0 * eos_ao3 * &
                (eosData(tempToUse+i)*eosData(tempToUse+i) / eosData(dens+i))**2
#else
           eosData(ded+i) = - eosData(eintRad+i) / eosData(dens+i) ! DEV: THIS dens?? - KW
#endif
        else                    !wants ions, I guess!
           eosData(ded+i) = 0.0
        end if
     end if
     if(maskPtr(EOS_DPT).or. .TRUE.) then
        if (wantComb) then
!!$           eosData(dpt+i) = eos_gasConstant*eosData(dens+i) * (eosData(zbar+i)+1.0) / eosData(abar+i)
           eosData(dpt+i) = eos_gasConstant*eosData(dens+i) / eosData(abar+i) &
                * (eosData(zbar+i) + 1.0 + &
                   eosData(tempToUse+i)*tabData(EOS_TABINT_DERIV_DT,EOS_TABVT_ZF,EOS_TAB_FOR_ELE))
#ifdef RAD_HERE
           eosData(dpt+i) = eosData(dpt+i) + 4.0 * eos_ao3 * eosData(tempToUse+i)**3
#else
           eosData(dpt+i) = eosData(dpt+i) + 4.0/3.0 * eosData(eintRad+i) * eosData(dens+i) / eosData(tempToUse+i)
#endif
        else if (wantEle) then
           eosData(dpt+i) = eos_gasConstant*eosData(dens+i) / eosData(abar+i) &
                * (eosData(zbar+i) + &
                   eosData(tempToUse+i)*tabData(EOS_TABINT_DERIV_DT,EOS_TABVT_ZF,EOS_TAB_FOR_ELE))
        else if (wantRad) then
#ifdef RAD_HERE
           eosData(dpt+i) = 4.0 * eos_ao3 * eosData(tempToUse+i)**3
#else
           eosData(dpt+i) = 4.0/3.0 * eosData(eintRad+i) * eosData(dens+i) / eosData(tempToUse+i)
#endif
        else                    !wants ions, I guess!
           eosData(dpt+i) = eos_gasConstant*eosData(dens+i) / eosData(abar+i)
        end if
     end if
     if(maskPtr(EOS_DPD).or. .TRUE.) then
        if (wantComb) then
!!$           eosData(dpd+i) = eos_gasConstant*eosData(tempToUse+i) * (eosData(zbar+i)+1.0) / eosData(abar+i)
           eosData(dpd+i) = eos_gasConstant*eosData(tempToUse+i) / eosData(abar+i) &
                * (eosData(zbar+i)+1+ eosData(dens+i)*tabData(EOS_TABINT_DERIV_DD,EOS_TABVT_ZF,EOS_TAB_FOR_ELE))
        else if (wantEle) then
           eosData(dpd+i) = eos_gasConstant*eosData(tempToUse+i) / eosData(abar+i) &
                * (eosData(zbar+i) + eosData(dens+i)*tabData(EOS_TABINT_DERIV_DD,EOS_TABVT_ZF,EOS_TAB_FOR_ELE))
        else if (wantRad) then
           eosData(dpd+i) = 0.0
        else                    !wants ions, I guess!
           eosData(dpd+i) = eos_gasConstant*eosData(tempToUse+i) / eosData(abar+i)
        end if
     end if
     if (.TRUE.) then
        if (wantComb) then
        else if (wantEle) then
           eosData(gamc+i) = eosData(tempToUse+i)/eosData(dens+i) * &
                             eosData(dpt+i)**2 /  eosData(det+i) &
                           + eosData(dens+i) * eosData(dpd+i)
        else if (wantRad) then
           eosData(gamc+i) = eosData(tempToUse+i)/eosData(dens+i) * &
                             eosData(dpt+i)**2 /  eosData(det+i) &
                           + eosData(dens+i) * eosData(dpd+i)
        else                    !wants ions, I guess!
           eosData(gamc+i) = eosData(tempToUse+i)/eosData(dens+i) * &
                             eosData(dpt+i)**2 /  eosData(det+i) &
                           + eosData(dens+i) * eosData(dpd+i)
        end if
     end if
  end do

  ! density, temperature taken as input
  select case (mode)
  case (MODE_DENS_TEMP)         !always wantComb...
     if (wantComb) &
          eosData(pres+ilo:pres+ihi) = eos_gasConstant*eosData(dens+ilo:dens+ihi) * &
            eosData(temp+ilo:temp+ihi) / eosData(abar+ilo:abar+ihi) * Zp
!!!$!!!     eosData(gamc+ilo:gamc+ihi) = 5./3.

     eosData(entr+ilo:entr+ihi) = 0.0

  case (MODE_DENS_TEMP_ION)
     if(maskPtr(EOS_PRESION)) &
          eosData(presIon+ilo:presIon+ihi) = eos_gasConstant * &
            eosData(dens+ilo:dens+ihi) * &
            eosData(tempIon+ilo:tempIon+ihi) / eosData(abar+ilo:abar+ihi)

     if(.NOT.maskPtr(EOS_PRESION))call Driver_abortFlash('eos_tabIonmix needs EOS_PRESION in mask!') 
     eosData(gamc+ilo:gamc+ihi) = eosData(gamc+ilo:gamc+ihi) / eosData(presIon+ilo:presIon+ihi)

     if(maskPtr(EOS_EINTION)) &
          eosData(eintIon+ilo:eintIon+ihi) = ggprodIon(ilo:ihi) * eosData(tempIon+ilo:tempIon+ihi) &
          / eosData(abar+ilo:abar+ihi)

     eosData(entr+ilo:entr+ihi) = 0.0
#ifdef DEBUG_EOS
     print*,'presIon by tabIonmix: ', eosData(presIon+ilo:presIon+ihi)
#endif


  case (MODE_DENS_TEMP_ELE)
     if(maskPtr(EOS_PRESELE)) &
          eosData(presEle+ilo:presEle+ihi) = eos_gasConstant * &
          eosData(dens+ilo:dens+ihi) * &
          eosData(tempEle+ilo:tempEle+ihi) / eosData(abar+ilo:abar+ihi) * dynamicZ

     if(.NOT.maskPtr(EOS_PRESELE))call Driver_abortFlash('eos_tabIonmix needs EOS_PRESELE in mask!') 
     eosData(gamc+ilo:gamc+ihi) = eosData(gamc+ilo:gamc+ihi) / eosData(presEle+ilo:presEle+ihi)
!!!$!!!     eosData(gamc+ilo:gamc+ihi) = 5./3.

!!$     if(maskPtr(EOS_EINTELE)) &
!!$          eosData(eintEle+ilo:eintEle+ihi) = ggprodEle * eosData(tempEle+ilo:tempEle+ihi) &
!!$          / eosData(abar+ilo:abar+ihi) * dynamicZ

     eosData(entr+ilo:entr+ihi) = 0.0
!!$     if(maskPtr(EOS_ENTRELE)) &
!!$          call setEntropy(entrele,tempEle)
     eosData(entrEle+ilo:entrEle+ihi) = 0.0
#ifdef DEBUG_EOS
     print*,'presEle by tabIonmix: ', eosData(presEle+ilo:presEle+ihi)
#endif


#if 0
  case (MODE_DENS_TEMP_RAD)
     if(maskPtr(EOS_PRESRAD)) &
          eosData(presRad+ilo:presRad+ihi)    = eos_ao3 * eosData(tempRad+ilo:tempRad+ihi)**4
     if(.NOT.maskPtr(EOS_PRESRAD))call Driver_abortFlash('eos_tabIonmix needs EOS_PRESRAD in mask!') 
     eosData(gamc+ilo:gamc+ihi) = eosData(gamc+ilo:gamc+ihi) / eosData(presRad+ilo:presRad+ihi)
     if(maskPtr(EOS_EINTRAD)) &
          eosData(eintRad+ilo:eintRad+ihi)    = 3.0e0 * eosData(presRad+ilo:presRad+ihi) / eosData(dens+ilo:dens+ihi)
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

     eosData(eintIon+ilo:eintIon+ihi) = ggprodIon(ilo:ihi) * eosData(tempIon+ilo:tempIon+ihi) &
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
        call eos_newtonRaphson(vecLen, MODE_DENS_EI, .TRUE., maskPtr,cMask=overallCMask)

        ! Next lines moved from end of eos_newtonRaphson; pres changed - KW 2010-12-13
        eosData(pres+ilo:pres+ihi) = ( &
             eosData(presIon+ilo:presIon+ihi)+ &
             eosData(presEle+ilo:presEle+ihi)+ &
             eosData(presRad+ilo:presRad+ihi))
#ifdef EOS_FORCE_2T
        eosData(gamc+ilo:gamc+ihi)=eos_gammaEle
#else
        call setCombinedGamc(gamc,coupled=.FALSE.)
!        eosData(gamc+ilo:gamc+ihi)=gamcRow(1:rowLen) !replace with expression based on component states!
#endif
        eosData(entr+ilo:entr+ihi)=stotRow(1:rowLen) !(replace with expression based on component states)

     end if

    !! Note that we make here a choice for the arbitrary additive constant in the electron entropy 
     if(maskPtr(EOS_ENTRELE)) &
          call setEntropy(entrele,tempEle)
#endif
     
        
  case (MODE_DENS_TEMP_EQUI,MODE_DENS_TEMP_MAT_EQUI)
     !Equalize temperature.  Each component uses same temperature.
     eosData(tempIon+ilo:tempIon+ihi) = eosData(temp+ilo:temp+ihi)
     eosData(tempEle+ilo:tempEle+ihi) = eosData(temp+ilo:temp+ihi)
     eosData(tempRad+ilo:tempRad+ihi) = eosData(temp+ilo:temp+ihi)

#ifdef RAD_HERE
     eosData(pres+ilo:pres+ihi) = eos_gasConstant*eosData(dens+ilo:dens+ihi) * &
          eosData(temp+ilo:temp+ihi) / eosData(abar+ilo:abar+ihi) * Zp         &
          +  eos_ao3 * eosData(tempRad+ilo:tempRad+ihi)**4
#else
     eosData(pres+ilo:pres+ihi) = eos_gasConstant*eosData(dens+ilo:dens+ihi) * &
          eosData(temp+ilo:temp+ihi) / eosData(abar+ilo:abar+ihi) * Zp         &
          +  eosData(presRad+ilo:presRad+ihi)
#endif
     eosData(presIon+ilo:presIon+ihi) = eos_gasConstant * &
          eosData(dens+ilo:dens+ihi) * &
          eosData(temp+ilo:temp+ihi) / eosData(abar+ilo:abar+ihi)
     eosData(presEle+ilo:presEle+ihi) = eos_gasConstant * &
          eosData(dens+ilo:dens+ihi) * &
          eosData(temp+ilo:temp+ihi) / eosData(abar+ilo:abar+ihi) * dynamicZ


     if (ANY(maskPtr((/EOS_EINTELE,EOS_EINTION,EOS_EINTRAD/)) .EQV. .FALSE.)) then
        call Driver_abortFlash("[eos_idealGamma] cannot calculate MODE_DENS_TEMP_EQUI without component eint masks.&
             & Set mask appropriately.")
     end if
     eosData(eintIon+ilo:eintIon+ihi) = ggprodIon(ilo:ihi) * eosData(temp+ilo:temp+ihi) &
          / eosData(abar+ilo:abar+ihi)
     eosData(eintEle+ilo:eintEle+ihi) =eosData(eint+ilo:eint+ihi) &
          - eosData(eintIon+ilo:eintIon+ihi) &
          - eosData(eintRad+ilo:eintRad+ihi)

#ifdef RAD_HERE
     if (maskPtr(EOS_PRESRAD)) &
          eosData(presRad+ilo:presRad+ihi)    = eos_ao3 * eosData(tempRad+ilo:tempRad+ihi)**4
#endif

     ! TESTING... :
     eosData(pres+ilo:pres+ihi) = ( &
          eosData(presIon+ilo:presIon+ihi)+ &
          eosData(presEle+ilo:presEle+ihi)+ &
          eosData(presRad+ilo:presRad+ihi))
     eosData(eint+ilo:eint+ihi) = ( &
          eosData(eintIon+ilo:eintIon+ihi)+ &
          eosData(eintEle+ilo:eintEle+ihi)+ &
          eosData(eintRad+ilo:eintRad+ihi))

     eosData(entr+ilo:entr+ihi) = 0.0

!!$     if(maskPtr(EOS_ENTRELE)) &
!!$          call setEntropy(entrele,temp)


9999 format(' E.o.MODE_DENS_TEMP_EQUI e:',1P,(7(1x,G20.6)))
     write(*,9999) (eosData(eint+i),eosData(eintIon+i),eosData(eintEle+i), eosData(eintRad+i), eosData(tempToUse+i),&
          eosData(dens+i),eosData(dens+i)/eosData(abar+i)*eos_avo ,i=ilo,ihi)
9998 format(' E.o.MODE_DENS_TEMP_EQUI p:',1P,4(1x,G20.6))
     write(*,9998) (eosData(pres+i),eosData(presIon+i),eosData(presEle+i), eosData(presRad+i),i=ilo,ihi)
#if(0)
!!$9997 format(' E.o.MODE_DENS_TEMP_EQUI Z:',1P,4(1x,G20.6))
!!$     write(*,9997) (eosData(zbar+i),i=ilo,ihi)
#endif
9996 format(' E.o.MODE_DENS_TEMP_EQUI det',1P,4(G20.6,1x))
     write(*,9996) (eosData(det+i),i=ilo,ihi)


#if 0
  case (MODE_DENS_TEMP_ALL)
     eosData(pres+ilo:pres+ihi) = eos_gasConstant*eosData(dens+ilo:dens+ihi) * &
          eosData(temp+ilo:temp+ihi) / eosData(abar+ilo:abar+ihi) * Zp

     eosData(presIon+ilo:presIon+ihi) = eos_gasConstant * &
          eosData(dens+ilo:dens+ihi) * &
          eosData(tempIon+ilo:tempIon+ihi) / eosData(abar+ilo:abar+ihi)
     eosData(presEle+ilo:presEle+ihi) = eos_gasConstant * &
          eosData(dens+ilo:dens+ihi) * &
          eosData(tempEle+ilo:tempEle+ihi) / eosData(abar+ilo:abar+ihi) * dynamicZ

     eosData(eint+ilo:eint+ihi) = ggprod(ilo:ihi) * eosData(temp+ilo:temp+ihi) &
          / eosData(abar+ilo:abar+ihi) * Zp
     eosData(eintIon+ilo:eintIon+ihi) = ggprodIon(ilo:ihi) * eosData(tempIon+ilo:tempIon+ihi) &
          / eosData(abar+ilo:abar+ihi)
     eosData(eintEle+ilo:eintEle+ihi) = ggprodEle * eosData(tempEle+ilo:tempEle+ihi) &
          / eosData(abar+ilo:abar+ihi) * dynamicZ

     eosData(presRad+ilo:presRad+ihi)    = eos_ao3 * eosData(tempRad+ilo:tempRad+ihi)**4

     eosData(eintRad+ilo:eintRad+ihi)    = 3.0e0 * eosData(presRad+ilo:presRad+ihi) / eosData(dens+ilo:dens+ihi)

     eosData(entr+ilo:entr+ihi) = 0.0

     if(maskPtr(EOS_ENTRELE)) &
          call setEntropy(entrele,tempEle)
#endif

  ! unrecognized or unsupported value for mode
  case default
     if (eos_meshMe==MASTER_PE) then
        print*,"eos_tabIonmix3T: invalid mode=",mode
        print*,"The Tabulated Eos implementation needs to be combined with something like the"
        print*,"Multitype Eos implementation in order to support modes that take energies or"
        print*,"pressures as input."
     end if
     call Driver_abortFlash("[eos_tabIonmix3T] Unrecognized input mode given to eos_tabIonmix")
  end select


  if(useNRForGamcCombined) then
     call Driver_abortFlash('Cannot set eosData(gamc+1...) from gamcRow w.d.n.e.!')
  else if (doCoupledForGamc .AND. (mode.NE.MODE_DENS_TEMP)) then
     call setCombinedGamc(gamc,coupled=.TRUE.)
  else if (doUncoupledForGamc) then
     call setCombinedGamc(gamc,coupled=.FALSE.)
  else
     ! leave as set above - KW
  end if


  if(present(mask)) then

     if(mask(EOS_DPT)) then
        dpt = (EOS_DPT-1)*vecLen
        if(useNRForAddtlCombined) then
           call Driver_abortFlash('Cannot set eosData(dpt+ilo...) from dptRow w.d.n.e.!')
!!           eosData(dpt+ilo:dpt+ihi) = dptRow(1:rowLen)
        else if (doUncoupledForAddtl) then
           eosData(dpt+ilo:dpt+ihi) = &
                (  eos_gasConstant * eosData(dens+ilo:dens+ihi) * &
                 ( eosData(tempIon+ilo:tempIon+ihi) + eosData(tempEle+ilo:tempEle+ihi)*dynamicZ)&
                  / eosData(abar+ilo:abar+ihi)  + &
                 4*eos_ao3 ) / eosData(tempEle+ilo:tempEle+ihi)
        else
           eosData(dpt+ilo:dpt+ihi) = eos_gasConstant*eosData(dens+ilo:dens+ihi) &
                * Zp / eosData(abar+ilo:abar+ihi)
        end if
     end if
     if(mask(EOS_DPD)) then
        dpd = (EOS_DPD-1)*vecLen
!!!!!!!!        eosData(dpd+ilo:dpd+ihi) = eos_gasConstant*eosData(temp+ilo:temp+ihi) / eosData(abar+ilo:abar+ihi)
     end if
     if(mask(EOS_DET))then
        det = (EOS_DET-1)*vecLen
        if(useNRForAddtlCombined) then
           call Driver_abortFlash('Cannot set eosData(det+ilo...) from detRow w.d.n.e.!')
!!           eosData(det+ilo:det+ihi) = detRow(1:rowLen)
#if(0)
        else if (.FALSE.) then
           call averageVectRByPressureFraction(det, &
                eosData(eintIon+ilo:eintIon+ihi) / eosData(tempIon+ilo:tempIon+ihi), &
                eosData(eintEle+ilo:eintEle+ihi) / eosData(tempEle+ilo:tempEle+ihi), &
                4 * eosData(eintRad+ilo:eintRad+ihi) / eosData(tempRad+ilo:tempRad+ihi))
#endif
        else if ((doCoupledForGamc .AND. (mode==MODE_DENS_TEMP_EQUI))) then
           ! NOP - already set above
        else if ((doCoupledForGamc .AND. (mode==MODE_DENS_TEMP_MAT_EQUI))) then
           ! NOP - already set above
        else if ((doCoupledForGamc .AND. (mode.NE.MODE_DENS_TEMP)).OR. doUncoupledForAddtl) then
           eosData(det+ilo:det+ihi) = &
                eosData(eintIon+ilo:eintIon+ihi) / eosData(tempIon+ilo:tempIon+ihi) + &
                eosData(eintEle+ilo:eintEle+ihi) / eosData(tempEle+ilo:tempEle+ihi) + &
                4 * eosData(eintRad+ilo:eintRad+ihi) / eosData(tempRad+ilo:tempRad+ihi) 
        else if (mode==MODE_DENS_TEMP_ION) then
           eosData(det+ilo:det+ihi) = &
                eosData(eintIon+ilo:eintIon+ihi) / eosData(tempIon+ilo:tempIon+ihi)
        else if (mode==MODE_DENS_TEMP_ELE) then
           ! Already done above from table lookup!
!!$           print*,'Duh!'
!!$           eosData(det+ilo:det+ihi) = &
!!$                eosData(eintEle+ilo:eintEle+ihi) / eosData(tempEle+ilo:tempEle+ihi)
        else if (mode==MODE_DENS_TEMP_RAD) then
           eosData(det+ilo:det+ihi) = &
                4 * eosData(eintRad+ilo:eintRad+ihi) / eosData(tempRad+ilo:tempRad+ihi)
!!$       else                     !laredy set above!
!!$
!!$           call Driver_abortFlash('Cannot set eosData(det+ilo...) from ggprod w.d.n.e.!')
!!$!!           eosData(det+ilo:det+ihi) = ggprod(ilo:ihi) * Zp / eosData(abar+ilo:abar+ihi)
        end if
     end if
!!$     if(mask(EOS_DED))then 
!!$        ded = (EOS_DED-1)*vecLen
!!$        eosData(ded+ilo:ded+ihi) = 0.
!!$     end if

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
              call Driver_abortFlash('Cannot set eosData(c_v+ilo...) from cvRow w.d.n.e.!')
!!              eosData(c_v+ilo:c_v+ihi) = cvRow(1:rowLen)
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
              call Driver_abortFlash('Cannot set eosData(c_p+ilo...) from cpRow w.d.n.e.!')
!!              eosData(c_p+ilo:c_p+ihi) = cpRow(1:rowLen)
           else if (doUncoupledForAddtl) then
              eosData(c_p+ilo:c_p+ihi) = eosData(game+ilo:game+ihi)*eosData(c_v+ilo:c_v+ihi)
           else
              eosData(c_p+ilo:c_p+ihi) = eosData(game+ilo:game+ihi)*eosData(c_v+ilo:c_v+ihi)
           end if
        else
           call Driver_abortFlash("[eos_idealGamma] cannot calculate C_P without C_V and DET.  Set mask appropriately.")
        end if
        
     end if

!!$     if(mask(EOS_CVION))then
!!$           c_v = (EOS_CVION-1)*vecLen
!!$           call Driver_abortFlash('Cannot set eosData(cv+ilo...) for EOS_CVION from ggprod w.d.n.e.!')
!!$!!           eosData(c_v+ilo:c_v+ihi) = ggprod(ilo:ihi)  / eosData(abar+ilo:abar+ihi)
!!$     end if
!!$     if(mask(EOS_CVELE))then
!!$           c_v = (EOS_CVELE-1)*vecLen
!!$           call Driver_abortFlash('Cannot set eosData(cv+ilo...) for EOS_CVELE from ggprod w.d.n.e.!')
!!$!!           eosData(c_v+ilo:c_v+ihi) = ggprod(ilo:ihi) * dynamicZ / eosData(abar+ilo:abar+ihi)
!!$     end if
  end if


  if (arrayBoundHackMode .eqv. .true.) then
     deallocate(maskPtr)
  end if

#ifdef DEBUG_EOS
  print*,'gamc a.e.o.tabIonmix: ', eosData(gamc+ilo:gamc+ihi)
#endif

  return

contains
  subroutine setCombinedGamc(gamcCombined,coupled)
    integer, intent(in) :: gamcCombined
    logical, intent(in) :: coupled

#if 0
    if (coupled) then
       eosData(gamcCombined+ilo:gamcCombined+ihi) = &
         ( gamM1Ion(ilo:ihi)*eosData(presIon+ilo:presIon+ihi)+ &
           eos_gammam1Ele*eosData(presEle+ilo:presEle+ihi)+ &
           eos_gammam1Rad*eosData(presRad+ilo:presRad+ihi) ) / eosData(pres+ilo:pres+ihi)
       eosData(gamcCombined+ilo:gamcCombined+ihi) = &
            1.0 + 1.0/eosData(gamcCombined+ilo:gamcCombined+ihi) 
    else
       eosData(gamcCombined+ilo:gamcCombined+ihi) = &
         ( gamIon(ilo:ihi)*eosData(presIon+ilo:presIon+ihi)+ &
           eos_gammaEle*eosData(presEle+ilo:presEle+ihi)+ &
           eos_gammaRad*eosData(presRad+ilo:presRad+ihi) ) / eosData(pres+ilo:pres+ihi)
    end if
#else
    eosData(gamcCombined+ilo:gamcCombined+ihi) = 5./3.
#endif

  end subroutine setCombinedGamc
end subroutine eos_tabIonmix
