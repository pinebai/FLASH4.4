!!****if* source/physics/Eos/EosMain/multiTemp/Multitype/eos_multiTypeByTemp
!!
!! NAME
!!
!!  eos_multiTypeByTemp
!!
!! SYNOPSIS
!!
!!  call      eos_multiTypeByTemp(integer(IN) :: mode,
!!                integer(IN)                 :: vecLen,
!!                real(INOUT)                 :: eosData(vecLen*EOS_NUM),
!!                integer(IN)                 :: vecBegin,
!!                integer(IN)                 :: vecEnd,
!!             optional,real(IN),dimension(*) :: massFrac(vecLen*NSPECIES),
!!             optional,target, logical(IN)   :: mask(EOS_VARS+1:EOS_NUM),
!!             optional, integer(IN)          :: componentMask(N_EOS_TEMP) )
!!
!! DESCRIPTION
!!
!!  Call Eos implementation routines of (possibly) different types and combine
!!  the results.
!!
!!  This subroutine is called from the Multitype implementation of Eos.F90,
!!  potentially from a Newton-Raphson loop for finding the right temperatures
!!  to satisfy energy requirements.
!!
!!  The model underlying THIS implementation (KW July 2011) is that various
!!  materials that co-occur in the same cells don't notice each other.
!!  Each material is handled as if it filled the cell (fully).
!!  Additive quantities are then added up to form the results returned.
!!
!!  This implementation is typically called for ions and electrons separately,
!!  and not necessarily at the same temperature.  While contributions from
!!  different species in the FLASH sence, i.e., different materials, are
!!  combined here, the caller is responsible for combining contributions
!!  of different components, i.e., ions and electrons (DEV: and radiation?).
!!
!! ARGUMENTS 
!! 
!!  mode :    Selects the mode of operation of the Eos unit.
!!            This routine will normally be called with a MODE_DENS_TEMP* mode.
!!            For current (FLASH4-beta) Multitype Eos.F90 purposes, the modes
!!            used will be MODE_DENS_TEMP_ELE and MODE_DENS_TEMP_ION if the mode
!!            for which Eos is called was MODE_DENS_TEMP_GATHER or
!!            MODE_DENS_EI_GATHER.
!!
!!  vecLen   : number of points (cells) for which the eosData array is sized.
!!             This is the maximum allowed value for vecEnd.
!!
!!  eosData  : This array is the data structure through which variable values are 
!!             passed in and out of the Eos routine. The arrays is sized as 
!!             EOS_NUM*vecLen. EOS_NUM, and individual input and output
!!             Eos variables are defined in Eos.h. The array is organizes such that
!!             the first 1:vecLen entries represent the first Eos variable, vecLen+1:
!!             2*vecLen represent the second Eos variable and so on.
!!             This is the same kind of packaging as for Eos.F90 etc.
!!
!!  vecBegin : Index of first cell in eosData to handle.
!!             Can be used to limit operation to a subrange of cells, untested.
!!             Must be greater than or equal to 1.
!!
!!  vecEnd   : Index of last cell in eosData to handle.
!!             Can be used to limit operation to a subrange of cells, untested.
!!             Must be less than or equal to vecLen.
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
!!  componentMask : A 3-element integer vector that contains 1 for componenets to
!!                  consider and 0 for components to ignore in this call, in the
!!                  order (ions, electrons, radiation).
!!                  This kind of duplicates information already present in the
!!                  mode argument (MODE_DENS_TEMP_ION vs. MODE_DENS_TEMP_ELE),
!!                  and may be useless and superfluous as well as redundant.
!!
!! NOTES
!!
!!  The global variable eos_type must be set to EOS_MTYPE in order for this routine to
!!  operate in Multitype mode; otherwise, it will just (pretty uselessly) pass the call
!!  on to the appropriate simpler Eos implementation.
!!
!!  NSPECIES is defined in Flash.h.
!!
!!  EOS_VARS and EOS_NUM  are defined in Eos.h.
!!  
!!  MODE_DENS_TEMP, MODE_DENS_EI, MODE_DENS_PRES, etc. are defined in constants.h.
!!
!!
!! SEE ALSO
!!
!!  Eos.h    defines the variables used.
!!  Eos_wrapped  sets up the data structure.
!!
!!***

subroutine eos_multiTypeByTemp(mode,vecLen,eosData,vecBegin,vecEnd,massFrac,mask,componentMask)

!==============================================================================
  use Multispecies_interface, ONLY : Multispecies_getProperty, &
       Multispecies_getPropertyVector, &
       Multispecies_getSumInv
  use Driver_interface, ONLY : Driver_abortFlash
  use Eos_data, ONLY : eos_meshMe, eos_type, eos_entrEleScaleChoice
  use eos_helmConstData, ONLY: eos_ao3, eos_avo, eos_kerg, eos_kergavo
!!$  use eos_vecData, ONLY:  detRow, dptRow !, dpdRow, dedRow, pelRow, neRow, etaRow, cvRow, cpRow
  use eos_localInterface, ONLY : eos_idealGamma3T, eos_mgamma, eos_helmholtz,&
       eos_tabulated, eos_tabIonmix
  implicit none
#include "constants.h"
#include "Flash.h"
#include "Eos.h"
#include "Multispecies.h"
  integer, INTENT(in) :: mode, vecLen
  real,INTENT(inout), dimension(EOS_NUM*vecLen) :: eosData 
  integer, INTENT(in) :: vecBegin,vecEnd
  real, optional, INTENT(in),dimension(:)    :: massFrac
  logical, optional, INTENT(in),target,dimension(EOS_VARS+1:EOS_NUM) :: mask
  integer, optional, INTENT(in),dimension(N_EOS_TEMP)::componentMask

  integer :: spec, specno
  integer :: ilo,ihi, rowLen
  integer :: passMode
  integer,dimension(N_EOS_TEMP) :: passCMask
  real, dimension(EOS_NUM*vecLen) :: curData, outData
  real, dimension(NSPECIES) :: Avect, Yvect
  real, dimension(NSPECIES,vecLen) :: Xmat, Ymat, Fmat !!!, Pmat
  real :: curA, curZ, curZMin
  real :: curGamma
  real :: abarInv
  integer :: specieStart, specieEnd
  integer :: dens, temp, pres, eint, abar, zbar, ekin
  integer :: gamc
  integer :: tempIon,tempEle,tempRad                                                
  integer :: eintIon,eintEle,eintRad, entrEle,entrRad, det, dst,cvion,cvele
  integer :: presIon,presEle,presRad, presComp, dpt
  integer :: i, v, offs
  integer :: subtype
  logical :: arrayBoundHackMode
  logical :: handleRadHere
  integer :: tempToUse
  integer :: tempToPrint
  integer :: speciesEosType, numberExistingSpecies

  logical,save,target, dimension(EOS_VARS+1:EOS_NUM) :: maskInternal
  data maskInternal / EOS_DERIVS * .FALSE. /
  logical,pointer, dimension(:) :: maskPtr

  logical :: presMFrac, presMask
  logical,parameter :: doRad = .TRUE.

  presMask = present(mask)
  presMFrac = present(massFrac)
  
  maskPtr => maskInternal
  if (presMask) then
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

  if (present(componentMask)) then
     passCMask = componentMask
  else
     passCMask = (/1,1,1/)
  end if
  if (passCMask(3) == 1 .AND. &
       (mode==MODE_DENS_TEMP_EQUI)) then !DEV: or... several others?
     passCMask(3) = 0
     handleRadHere = .TRUE.
  else
     handleRadHere = .FALSE.
  end if
     

!!$  if (present(vecBegin)) then
     ilo = vecBegin
!!$  else
!!$     ilo = 1
!!$  end if
!!$  if (present(vecEnd)) then
     ihi = vecEnd
!!$  else
!!$     ihi = vecLen
!!$  end if
!!$  rowLen = ihi - ilo + 1
!!$#ifdef DEBUG_EOS
!!$  if (ilo < 1 .OR. ilo > vecLen) then
!!$     print*,'[eos_multiTypeByTemp] ilo is',ilo
!!$     call Driver_abortFlash("[eos_multiTypeByTemp] invalid ilo")
!!$  end if
!!$  if (ihi < 1 .OR. ihi > vecLen) then
!!$     print*,'[eos_multiTypeByTemp] ihi is',ihi
!!$     call Driver_abortFlash("[eos_multiTypeByTemp] invalid ihi")
!!$  end if
!!$  if (rowLen < 0 .OR. rowLen > vecLen) then
!!$     print*,'[eos_multiTypeByTemp] rowLen is',rowLen
!!$     call Driver_abortFlash("[eos_multiTypeByTemp] invalid rowLen")
!!$  end if
!!$#endif
!!$  if (rowLen == 0) then
!!$     print*,'[eos_multiTypeByTemp] rowLen is 0.'
!!$  end if

  pres = (EOS_PRES-1)*vecLen
  dens = (EOS_DENS-1)*vecLen
  temp = (EOS_TEMP-1)*vecLen
  eint = (EOS_EINT-1)*vecLen
  gamc = (EOS_GAMC-1)*vecLen
  abar = (EOS_ABAR-1)*vecLen
  zbar = (EOS_ZBAR-1)*vecLen
  ekin = (EOS_EKIN-1)*vecLen
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
  det = (EOS_DET-1)*vecLen
  dpt = (EOS_DPT-1)*vecLen
  dst = (EOS_DST-1)*vecLen
  cvion = (EOS_CVION-1)*vecLen
  cvele = (EOS_CVELE-1)*vecLen

  select case (mode)
  ! density, temperature taken as input
  case (MODE_DENS_TEMP)
  end select


  if (presMask) then
     ! nothing special
  end if

  select case(eos_type)
  case(EOS_GAM)
     call eos_idealGamma3T(mode, vecLen, eosData,vecBegin,vecEnd, &
          eosType=eos_type, mask=mask, componentMask=componentMask)
  case(EOS_MGAM)
     call eos_mgamma(mode, vecLen, eosData, vecBegin,vecEnd, mask=mask)
  case(EOS_HLM)
     call eos_helmholtz(mode, vecLen, eosData, mask=mask)
  case(EOS_TAB)
     call eos_tabulated(mode, vecLen, eosData, mask=mask)
  case(EOS_MTYPE)               !do this, temporarily...
     if (.NOT. presMFrac) then
        call Driver_abortFlash('eos_multiTypeByTemp needs mass fractions!')
     end if
     outData(:) = 0.0
     curData(:) = eosData(:)
     
     do i = vecBegin,vecEnd
        specieStart = (i-1)*NSPECIES + 1
        specieEnd = i*NSPECIES
        call Multispecies_getSumInv(A, abarInv ,massFrac(specieStart:specieEnd))
        eosData(abar+i) = 1.0/abarInv
     end do

     call Multispecies_getPropertyVector(A,Avect)
     do specno = 1,NSPECIES
        Xmat(specno,vecBegin:vecEnd) = massFrac((vecBegin-1)*NSPECIES+specno:(vecEnd-1)*NSPECIES+specno:NSPECIES)
        Ymat(specno,vecBegin:vecEnd) = Xmat(specno,vecBegin:vecEnd) / Avect(specno)
        Fmat(specno,vecBegin:vecEnd) = Ymat(specno,vecBegin:vecEnd) * eosData(abar+vecBegin:abar+vecEnd)
!!!        Pmat(specno,vecBegin:vecEnd) = 0.0
     end do

!!$     print*,'Ymat:',Ymat

     numberExistingSpecies = 0
     do spec=SPECIES_BEGIN,SPECIES_END
        specno = spec - SPECIES_BEGIN + 1
        if (ALL(massFrac((vecBegin-1)*NSPECIES+specno:(vecEnd-1)*NSPECIES+specno:NSPECIES).LE. 1.0e-20)) then
#ifdef DEBUG_EOS
           print*,' skipping species #',specno,' for',&
                size(massFrac((vecBegin-1)*NSPECIES+specno:(vecEnd-1)*NSPECIES+specno:NSPECIES)),&
                ' cells since all massfracs are (nearly) 0.'
#endif
        else
!!$           print*,'NOT Skipping species #',specno,' for',&
!!$                size(massFrac((vecBegin-1)*NSPECIES+specno:(vecEnd-1)*NSPECIES+specno:NSPECIES)),&
!!$                ' cells:',&
!!$                massFrac((vecBegin-1)*NSPECIES+specno:(vecEnd-1)*NSPECIES+specno:NSPECIES)

           numberExistingSpecies = numberExistingSpecies + 1
           call Multispecies_getProperty(spec,MS_EOSTYPE,speciesEosType)
           do v=1,EOS_NUM
              offs=(v-1)*vecLen
              curData(offs+vecBegin:offs+vecEnd) = eosData(offs+vecBegin:offs+vecEnd)
           end do
           curData(dens+vecBegin:dens+vecEnd) = eosData(dens+vecBegin:dens+vecEnd) * Xmat(specno,vecBegin:vecEnd)

           select case (speciesEosType)
           case (EOS_GAM)
              call Multispecies_getProperty(spec,GAMMA,curGamma) !to delete!?
              call Multispecies_getProperty(spec,A    ,curA)
              call Multispecies_getProperty(spec,Z    ,curZ)
              call Multispecies_getProperty(spec,MS_ZMIN,curZMin)
              curZ = max(curZ,curZMin)
              curData(gamc+vecBegin:gamc+vecEnd) = curGamma !to delete!?
              curData(abar+vecBegin:abar+vecEnd) = curA
              curData(zbar+vecBegin:zbar+vecEnd) = curZ

              call eos_idealGamma3T(mode, vecLen, curData,vecBegin,vecEnd, &
                   eosType=EOS_GAM, material=spec, massFrac=massFrac, mask=mask, componentMask=passCMask)
           case (EOS_TAB)
              call Multispecies_getProperty(spec,A    ,curA)
              call Multispecies_getProperty(spec,MS_ZMIN,curZMin)
              call Multispecies_getProperty(spec,MS_EOSSUBTYPE,subtype)
              curData(abar+vecBegin:abar+vecEnd) = curA
              curData(zbar+vecBegin:zbar+vecEnd) = max(curZMin,0.0)
              if (doRad) then
                 tempToUse = temp !DEV: for now
                 do i = vecBegin,vecEnd
                    curData(eintRad+i) = 3.0e0 * eos_ao3 * curData(tempToUse+i)**4 / curData(dens+i) &
                         * Xmat(specno,i)
                    curData(presRad+i) = eos_ao3 * curData(tempToUse+i)**4  * Xmat(specno,i)
                 end do
              end if

              if (mode==MODE_DENS_TEMP_EQUI) then
                 passMode = MODE_DENS_TEMP_MAT_EQUI
              else
                 passMode = mode
              end if
              call eos_tabIonmix(passMode, vecLen, curData,vecBegin,vecEnd, &
                   eosType=EOS_TAB, subtype=subtype, material=spec, mask=mask)
           case default
              if (eos_meshMe==MASTER_PE) print*,'[eos_multiTypeByTemp] Invalid EOS type',speciesEosType, &
                      ' for species',spec
              call Driver_abortFlash('[eos_multiTypeByTemp] Invalid per-species EOS type')
           end select

           if (handleRadHere) then
              tempToUse = temp !DEV: for now
              do i = vecBegin,vecEnd
                 curData(eintRad+i) = 3.0e0 * eos_ao3 * curData(tempToUse+i)**4 / curData(dens+i) &
                      * Xmat(specno,i)
                 curData(presRad+i) = eos_ao3 * curData(tempToUse+i)**4  * Xmat(specno,i)
                 if (mask(EOS_DET)) &
                    curData(det+i)     = curData(det+i) + 12.0 * eos_ao3 * curData(tempToUse+i)**3 / curData(dens+i) &
                      * Xmat(specno,i)
                 if (mask(EOS_DPT)) &
                    curData(dpt+i)     = curData(dpt+i) + 4.0  * eos_ao3 * curData(tempToUse+i)**3 &
                      * Xmat(specno,i)
              end do
           end if
           call eos_rescaleCurData(curData,mode,maskPtr)
           call eos_accumData(outData,curData,mode)
        end if
     end do
     call eos_rescaleOutData(outData,mode)
     outData(abar+vecBegin:abar+vecEnd) = eosData(abar+vecBegin:abar+vecEnd)
     outData(ekin+vecBegin:ekin+vecEnd) = eosData(ekin+vecBegin:ekin+vecEnd)
     do v=1,EOS_NUM
        offs=(v-1)*vecLen
        if (v>EOS_VARS) then
           if (present(mask)) then
              if(.NOT. mask(v)) cycle
           else
              cycle
           end if
        end if
        eosData(offs+vecBegin:offs+vecEnd) = outData(offs+vecBegin:offs+vecEnd)
     end do
!!$     detRow(vecBegin:vecEnd) = eosData(det+vecBegin:det+vecEnd)
!!$     dptRow(vecBegin:vecEnd) = eosData(dpt+vecBegin:dpt+vecEnd)

  case default
     if (eos_meshMe==MASTER_PE) print*,'[eos_multiTypeByTemp] eos_type has invalid value',eos_type
     call Driver_abortFlash('[eos_multiTypeByTemp] invalid eos_type')
  end select

#ifdef DEBUG_EOS
9999 format(' E.o.multiTypeByTemp     e:',1P,(7(1x,G20.6)))
     write(*,9999) (eosData(eint+i),eosData(eintIon+i),eosData(eintEle+i), eosData(eintRad+i), eosData(temp+i),&
          eosData(dens+i),eosData(dens+i)/eosData(abar+i)*eos_avo ,i=ilo,ihi)
     select case (mode)
        case(MODE_DENS_TEMP_ION)
           tempToPrint = tempIon
        case(MODE_DENS_TEMP_ELE)
           tempToPrint = tempEle
        case(MODE_DENS_TEMP_RAD)
           tempToPrint = tempRad
        case default
           tempToPrint = temp
        end select
9998 format(' E.o.multiTypeByTemp     p:',1P,5(1x,G20.6))
     write(*,9998) (eosData(pres+i),eosData(presIon+i),eosData(presEle+i), eosData(presRad+i), eosData(tempToPrint+i),&
         i=ilo,ihi)
#if(0)
!!$9997 format(' E.o.multiTypeByTemp     Z:',1P,4(1x,G20.6))
!!$     write(*,9997) (eosData(zbar+i),i=ilo,ihi)
#endif
9996 format(' E.o.multiTypeByTemp     det',1P,4(G20.6,1x))
     write(*,9996) (eosData(det+i),i=ilo,ihi)
#endif


  if (arrayBoundHackMode .eqv. .true.) then
     deallocate(maskPtr)
  end if

  return

contains
  subroutine eos_rescaleCurData(curData, mode, mask)
    real,intent(INOUT),  dimension(EOS_NUM*vecLen) :: curData
    integer, intent(in) :: mode
    logical, INTENT(in),dimension(EOS_VARS+1:EOS_NUM) :: mask
    real,               dimension(vecLen) :: presRelevant
    real    :: prefactor(vecBegin:vecEnd)
    integer :: ic

!    curData(gamc+vecBegin:gamc+vecEnd) = Fmat(specno,vecBegin:vecEnd) / (curData(gamc+vecBegin:gamc+vecEnd)-1.0)
    if (mode==MODE_DENS_TEMP_ION.OR.mode==MODE_DENS_TEMP_ELE.OR.mode==MODE_DENS_TEMP_RAD.OR.&
         mode==MODE_DENS_TEMP_COMP) then
       presRelevant = 0.0
       do ic=0,N_EOS_TEMP-1
          if (passCMask(ic+1) .ne. 0) then
             presRelevant(vecBegin:vecEnd) = presRelevant(vecBegin:vecEnd) + &
                  curData(presIon+ic*vecLen+vecBegin:presIon+ic*vecLen+vecEnd)
          end if

       end do
    else
       presRelevant(vecBegin:vecEnd) = curData(pres+vecBegin:pres+vecEnd)
    end if
#ifdef GAMC_AVERAGE_COUPLED
    curData(gamc+vecBegin:gamc+vecEnd) = &
         presRelevant(vecBegin:vecEnd) / (curData(gamc+vecBegin:gamc+vecEnd)-1.0)
#else
    curData(gamc+vecBegin:gamc+vecEnd) = &
         presRelevant(vecBegin:vecEnd) * curData(gamc+vecBegin:gamc+vecEnd)
#endif

!!$    curData(abar+vecBegin:abar+vecEnd) = curData(abar+vecBegin:abar+vecEnd) * Xmat(specno,vecBegin:vecEnd)
    curData(zbar+vecBegin:zbar+vecEnd) = curData(zbar+vecBegin:zbar+vecEnd) * Ymat(specno,vecBegin:vecEnd)

    curData(eint+vecBegin:eint+vecEnd) = curData(eint+vecBegin:eint+vecEnd) * Xmat(specno,vecBegin:vecEnd)
    if (mask(EOS_EINTION)) &
         curData(eintIon+vecBegin:eintIon+vecEnd) = curData(eintIon+vecBegin:eintIon+vecEnd) * Xmat(specno,vecBegin:vecEnd)
    if (mask(EOS_EINTELE)) &
         curData(eintEle+vecBegin:eintEle+vecEnd) = curData(eintEle+vecBegin:eintEle+vecEnd) * Xmat(specno,vecBegin:vecEnd)
    if (mask(EOS_EINTRAD)) &
         curData(eintRad+vecBegin:eintRad+vecEnd) = curData(eintRad+vecBegin:eintRad+vecEnd) * Xmat(specno,vecBegin:vecEnd)
    curData(entrEle+vecBegin:entrEle+vecEnd) = curData(entrEle+vecBegin:entrEle+vecEnd) * Xmat(specno,vecBegin:vecEnd)
!!$    print*,'EntrA:', curData(entrEle+vecBegin:entrEle+vecEnd)
    if (mask(EOS_ENTRRAD)) &
         curData(entrRad+vecBegin:entrRad+vecEnd) = curData(entrRad+vecBegin:entrRad+vecEnd) * Xmat(specno,vecBegin:vecEnd)
    if (mode==MODE_DENS_TEMP_ELE) then
       select case (eos_entrEleScaleChoice)
       case(1)
          prefactor = eos_kerg
       case(2,3)
          prefactor = eosData(zbar+vecBegin:zbar+vecEnd) / eosData(abar+vecBegin:abar+vecEnd) * eos_kergavo
       case(4,6)
          prefactor = eosData(zbar+vecBegin:zbar+vecEnd) / eosData(abar+vecBegin:abar+vecEnd) ! Ye
       case(5,7)
          prefactor = 1.0
       case default
          call Driver_abortFlash('eos_multiTypeByTemp - unsupported eos_entrEleScaleChoice value')
       end select
       select case (eos_entrEleScaleChoice)
       case(4,6)
          curData(entrEle+vecBegin:entrEle+vecEnd) = curData(entrEle+vecBegin:entrEle+vecEnd) + &
            curData(zbar+vecBegin:zbar+vecEnd) * log(Xmat(specno,vecBegin:vecEnd)) ! includes * Ymat(specno,vecBegin:vecEnd)
       case default
          curData(entrEle+vecBegin:entrEle+vecEnd) = curData(entrEle+vecBegin:entrEle+vecEnd) + &
            prefactor * Fmat(specno,vecBegin:vecEnd) * log(Xmat(specno,vecBegin:vecEnd))
       end select
!!$       print*,'EntrB:', curData(entrEle+vecBegin:entrEle+vecEnd)
    end if

    if (mask(EOS_DET)) &
         curData(det+vecBegin:det+vecEnd) = curData(det+vecBegin:det+vecEnd) * Xmat(specno,vecBegin:vecEnd)
    if (mask(EOS_DST)) &
         curData(dst+vecBegin:dst+vecEnd) = curData(dst+vecBegin:dst+vecEnd) * Xmat(specno,vecBegin:vecEnd)
    if (mask(EOS_CVION)) &
         curData(cvion+vecBegin:cvion+vecEnd) = curData(cvion+vecBegin:cvion+vecEnd) * Xmat(specno,vecBegin:vecEnd)
    if (mask(EOS_CVELE)) &
         curData(cvele+vecBegin:cvele+vecEnd) = curData(cvele+vecBegin:cvele+vecEnd) * Xmat(specno,vecBegin:vecEnd)
  end subroutine eos_rescaleCurData

  subroutine eos_accumData(cumData, curData, mode)
    real,intent(INOUT), dimension(EOS_NUM*vecLen) :: cumData
    real,intent(in),  dimension(EOS_NUM*vecLen) :: curData
    integer, intent(in) :: mode

    integer :: vv,offs

    do vv=1,EOS_NUM
       offs=(vv-1)*vecLen
       if (vv>EOS_VARS) then
          if (present(mask)) then
             if(.NOT. mask(vv)) cycle
          else
             cycle
          end if
       end if
       cumData(offs+vecBegin:offs+vecEnd) = cumData(offs+vecBegin:offs+vecEnd) + curData(offs+vecBegin:offs+vecEnd)
    end do
    cumData(temp+vecBegin:temp+vecEnd) = curData(temp+vecBegin:temp+vecEnd)
    cumData(tempIon+vecBegin:tempIon+vecEnd) = curData(tempIon+vecBegin:tempIon+vecEnd)
    cumData(tempEle+vecBegin:tempEle+vecEnd) = curData(tempEle+vecBegin:tempEle+vecEnd)
    cumData(tempRad+vecBegin:tempRad+vecEnd) = curData(tempRad+vecBegin:tempRad+vecEnd)
  end subroutine eos_accumData

  subroutine eos_rescaleOutData(outData, mode)
    real,intent(INOUT),  dimension(EOS_NUM*vecLen) :: outData
    integer, intent(in) :: mode
    real,               dimension(vecLen) :: presRelevant
    integer :: ic

    if (mode==MODE_DENS_TEMP_ION.OR.mode==MODE_DENS_TEMP_ELE.OR.mode==MODE_DENS_TEMP_RAD.OR.&
         mode==MODE_DENS_TEMP_COMP) then
       presRelevant = 0.0
       do ic=0,N_EOS_TEMP-1
          if(passCMask(ic+1) .ne. 0) then
             presRelevant(vecBegin:vecEnd) = presRelevant(vecBegin:vecEnd) + &
                  outData(presIon+ic*vecLen+vecBegin:presIon+ic*vecLen+vecEnd)
          end if
       end do
    else
       presRelevant(vecBegin:vecEnd) = outData(pres+vecBegin:pres+vecEnd)
    end if

!    outData(gamc+vecBegin:gamc+vecEnd) = 1.0 +  1.0/outData(gamc+vecBegin:gamc+vecEnd)
#ifdef GAMC_AVERAGE_COUPLED
    outData(gamc+vecBegin:gamc+vecEnd) = 1.0 +  presRelevant(vecBegin:vecEnd)/outData(gamc+vecBegin:gamc+vecEnd)
#else
    outData(gamc+vecBegin:gamc+vecEnd) = outData(gamc+vecBegin:gamc+vecEnd) / presRelevant(vecBegin:vecEnd)
#endif

    outData(zbar+vecBegin:zbar+vecEnd) = outData(zbar+vecBegin:zbar+vecEnd) * eosData(abar+vecBegin:abar+vecEnd)

  end subroutine eos_rescaleOutData

end subroutine eos_multiTypeByTemp
