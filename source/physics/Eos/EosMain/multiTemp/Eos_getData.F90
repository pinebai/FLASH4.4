!!****if* source/physics/Eos/EosMain/multiTemp/Eos_getData
!! NAME
!!
!!  Eos_getData
!! 
!! SYNOPSIS
!!
!!  call Eos_getData(  integer(IN) :: axis,
!!                     integer(IN) :: pos(MDIM),
!!                     integer(IN) :: vecLen,
!!                  real, pointer  :: solnData(:,:,:,:),
!!                     integer(IN) :: gridDataStruct,
!!                     real(OUT)   :: eosData(:),
!!            optional,real(OUT)   :: massFrac(:),
!!         optional,logical(INOUT) :: eosMask(EOS_VARS+1:) )
!!
!!
!!
!! DESCRIPTION
!!
!! Eos_getData gets data from a Grid data structure into an eosData array, for
!! passing to a subsequent Eos call. Additionally, mass fractions from the Grid
!! data structure may be extracted to the optional massFrac argument array.
!!
!! The Eos_wrapped function is provided for the user's convenience and acts as a simple
!! wrapper to the Eos interface. The Eos interface uses a single, flexible data
!! structure "eosData" to pass the thermodynamic quantities in and out of the
!! function (see Eos). The wrapper hides formation and use of eosData
!! from the users. The wrapper function uses the Eos_getData function to construct the 
!! data structure eosData.
!!
!! While Eos does not know anything about blocks, Eos_getData takes its
!! input thermodynamic state variables from a given block's storage area,
!! a vector at a time. It works by taking a selected vector of a block
!! described by the arguments axis, pos and vecLen.
!!
!!
!!  ARGUMENTS 
!!
!!   
!!   axis : the dimension of the vector in the block's storage
!!   pos  : the starting indices of the vector in the block. Note that the
!!          vector has to provide the starting indices for all dimensions
!!   vecLen : the length of the vector
!!   solnData: data from the current block; unmodified on return.
!!              various components (variables) of solnData will determine the
!!              contents of eosData.
!!   gridDataStruct : the relevant grid data structure, on whose data Eos was applied.
!!                    One of CENTER, FACEVAR{X,Y,Z}, GRIDVAR, defined in constants.h .
!!   eosData : the data structure native to the Eos unit; input and to Eos() as
!!             well as output from Eos().
!!   massFrac : this is an optional argument which is used when there is more 
!!              than one species in the simulation
!!   eosMask : if the caller passes in an eosMask array, Eos_getData may modify
!!             this mask for the following Eos() calls, probably removing
!!             some of the requested derived quantities.
!!
!!  EXAMPLE 
!!      if axis = IAXIS, pos(IAXIS)=1,pos(JAXIS)=1,pos(KAXIS)=1 and vecLen=4
!!      then Eos is to be applied to four cells in the first row along IAXIS
!!      of the lower left hand corner of the guard cells in the block. 
!!
!!      However if the value were
!!         pos(IAXIS)=iguard+1,
!!         pos(JAXIS)=jguard+1,
!!         pos(KAXIS)=kguard+1, vecLen = NYB, and axis = JAXIS
!!      then Eos is applied to the first column along Y axis in the interior of the block.
!!
!!  NOTES
!!
!!      This interface is called from Eos_wrappped, and is normally not called
!!      by user code directly.
!!
!!      The actual arguments in a call should match those used in a subsequent
!!      Eos_putData call that will extract results from the Eos() call out of
!!      the eosData array.
!!
!!      This interface is defined in Fortran Module 
!!      Eos_interface. All functions calling this routine should include
!!      a statement like
!!      use Eos_interface, ONLY : Eos_putData
!!
!!      This routine cannot use "INTERIOR" mode of indexing the range.  In the
!!      second example given above, although only the interior cells are being
!!      calculated with EOS, the range indices still must include the guard cells.
!!      See, for example, IsentropicVortex/Simulation_initBlock where the data is
!!      generated on INTERIOR cells with Grid_putRowData, but the same indices can't
!!      be used for the EOS call.
!!
!!  SEE ALSO
!!
!!     Eos
!!     Eos_putData
!!     Eos.h
!!
!!
!!***

! solnData depends on the ordering on unk
!!REORDER(4): solnData


subroutine Eos_getData(axis,pos,vecLen,solnData,gridDataStruct,eosData,massFrac, eosMask)

  use Driver_interface, ONLY: Driver_abortFlash
  use Eos_data, ONLY: eos_eintSwitch, eos_smalle, eos_smallEion, eos_smallEele, eos_smallErad, &
                      eos_smallT, eos_mapLookup

  implicit none
  
#include "Eos.h"
#include "Eos_map.h"
#include "constants.h"
#include "Flash.h"
  
  integer, intent(in) :: axis, vecLen, gridDataStruct
  integer, dimension(MDIM), intent(in) :: pos
  real, dimension(:),intent(OUT) :: eosData
  real,dimension(:),optional,intent(OUT) :: massFrac
  logical, optional, INTENT(INOUT),dimension(EOS_VARS+1:) :: eosMask     
  real, pointer:: solnData(:,:,:,:)


  integer :: i,j,k,n,m,pres,dens,gamc,temp,abar,zbar,eint,ekin,entr
  integer :: tempIon,tempEle,tempRad
  integer :: eintIon,eintEle,eintRad
  integer :: presIon,presEle,presRad
  integer :: pres_map,dens_map,gamc_map,temp_map,entr_map
  integer :: eint_map,ener_map, velx_map, vely_map, velz_map, sumy_map, ye_map
  integer :: temp1_map,temp2_map,temp3_map
  integer :: eint1_map,eint2_map,eint3_map
  integer :: pres1_map,pres2_map,pres3_map
  integer :: entrEle, sele_map
  integer :: entrRad, srad_map
  integer :: ib,ie,jb,je,kb,ke
  real :: kineticEnergy, internalEnergy
  real :: internalIonEnergy
  real :: internalEleEnergy
  real :: internalRadEnergy
!! ---------------------------------------------------------------------------------
  ! Test calling arguments

  ! Initializations:   grab the solution data from UNK and determine
  !   the length of the data being operated upon
  
  ib=pos(IAXIS)
  jb=pos(JAXIS)
  kb=pos(KAXIS)
  ie=pos(IAXIS)
  je=pos(JAXIS)
  ke=pos(KAXIS)
  select case(axis)
  case(IAXIS)
     ie=ie+vecLen-1
  case(JAXIS)
     je=je+vecLen-1
  case(KAXIS)
     ke=ke+vecLen-1
  end select
  ! These integers are indexes into the location in eosData just before the storage area for the appropriate variable.
  pres = (EOS_PRES-1)*vecLen
  dens = (EOS_DENS-1)*vecLen
  temp = (EOS_TEMP-1)*vecLen
  gamc = (EOS_GAMC-1)*vecLen
  eint = (EOS_EINT-1)*vecLen
  ekin = (EOS_EKIN-1)*vecLen
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
#ifdef EOS_ENTRELE
  entrEle = (EOS_ENTRELE-1)*vecLen
#endif
#ifdef EOS_ENTRRAD
  entrRad = (EOS_ENTRRAD-1)*vecLen
#endif

  pres_map = eos_mapLookup(EOSMAP_PRES,EOS_IN,gridDataStruct)
  dens_map = eos_mapLookup(EOSMAP_DENS,EOS_IN,gridDataStruct)
  temp_map = eos_mapLookup(EOSMAP_TEMP,EOS_IN,gridDataStruct)
  gamc_map = eos_mapLookup(EOSMAP_GAMC,EOS_IN,gridDataStruct)
  eint_map = eos_mapLookup(EOSMAP_EINT,EOS_IN,gridDataStruct)
  ener_map = eos_mapLookup(EOSMAP_ENER,EOS_IN,gridDataStruct)
  velx_map = eos_mapLookup(EOSMAP_VELX,EOS_IN,gridDataStruct)
  vely_map = eos_mapLookup(EOSMAP_VELY,EOS_IN,gridDataStruct)
  velz_map = eos_mapLookup(EOSMAP_VELZ,EOS_IN,gridDataStruct)
  sumy_map = eos_mapLookup(EOSMAP_SUMY,EOS_IN,gridDataStruct)
  ye_map = eos_mapLookup(EOSMAP_YE,EOS_IN,gridDataStruct)
  entr_map = eos_mapLookup(EOSMAP_ENTR,EOS_IN,gridDataStruct)
  temp1_map = eos_mapLookup(EOSMAP_TEMP1,EOS_IN,gridDataStruct)
  temp2_map = eos_mapLookup(EOSMAP_TEMP2,EOS_IN,gridDataStruct)
  temp3_map = eos_mapLookup(EOSMAP_TEMP3,EOS_IN,gridDataStruct)
  eint1_map = eos_mapLookup(EOSMAP_EINT1,EOS_IN,gridDataStruct)
  eint2_map = eos_mapLookup(EOSMAP_EINT2,EOS_IN,gridDataStruct)
  eint3_map = eos_mapLookup(EOSMAP_EINT3,EOS_IN,gridDataStruct)
  pres1_map = eos_mapLookup(EOSMAP_PRES1,EOS_IN,gridDataStruct)
  pres2_map = eos_mapLookup(EOSMAP_PRES2,EOS_IN,gridDataStruct)
  pres3_map = eos_mapLookup(EOSMAP_PRES3,EOS_IN,gridDataStruct)
  if (present(eosMask)) then
#ifdef EOS_ENTRELE
     sele_map = eos_mapLookup(EOSMAP_SELE,EOS_OUT,gridDataStruct)
     if (sele_map==NONEXISTENT) eosMask(EOS_ENTRELE) = .FALSE.
#endif
#ifdef EOS_ENTRRAD
     srad_map = eos_mapLookup(EOSMAP_SRAD,EOS_OUT,gridDataStruct)
     if (srad_map==NONEXISTENT) eosMask(EOS_ENTRRAD) = .FALSE.
#endif
  end if
  sele_map = eos_mapLookup(EOSMAP_SELE,EOS_IN,gridDataStruct)
  srad_map = eos_mapLookup(EOSMAP_SRAD,EOS_IN,gridDataStruct)
     

  if(present(massFrac)) then
     if(gridDataStruct /= SCRATCH) then
        m=1
        do k = kb,ke
           do j = jb,je
              do i = ib,ie
                 do n = SPECIES_BEGIN,SPECIES_END
                    massFrac(m) = solnData(n,i,j,k)
                    m=m+1
                 end do
              end do
           end do
        end do
     else
        m=1
        do k = kb,ke
           do j = jb,je
              do i = ib,ie
                 do n = SPECIES_BEGIN,SPECIES_END
                    massFrac(m) = solnData(i,j,k,n)
                    m=m+1
                 end do
              end do
           end do
        end do
     end if
  end if

  if(gridDataStruct /= SCRATCH) then
     n = 0
     do k = kb,ke
        do j = jb,je
           do i = ib,ie
              if (velx_map > 0 .AND. vely_map > 0 .AND. velz_map > 0) then
                 kineticEnergy  = 0.5*(solnData(velx_map,i,j,k)**2 + &
                                       solnData(vely_map,i,j,k)**2 + &
                                       solnData(velz_map,i,j,k)**2)
              else
                 kineticEnergy = 0.0
              end if

              n=n+1
              eosData(ekin+n) = kineticEnergy
        !! kineticEnergy holds velocity vector information -- 1/2 * Vmag**2
        !! internalEnergy holds eint (directly)  or energyTotal - ekinetic (calculated),
        !!          depending upon eintSwitch
              if(eint_map /= NONEXISTENT) then
                 internalEnergy  = solnData(eint_map,i,j,k)
                 if(ener_map /= NONEXISTENT) then
                 if ( solnData(ener_map,i,j,k) - kineticEnergy > max(eos_smalle, eos_eintSwitch*kineticEnergy)) then
                       internalEnergy = solnData(ener_map,i,j,k) - kineticEnergy
#ifdef DEBUG_EOS
                    else if (kineticEnergy .NE. 0.0) then
                       print*,'Eos_getD: Using EINT_VAR!',solnData(ener_map,i,j,k)-internalEnergy,kineticEnergy
#endif
                    end if
                 end if
              else
                 internalEnergy = solnData(ener_map,i,j,k)-kineticEnergy
              endif
#ifdef DEBUG_EOS
              if (internalEnergy<eos_smalle)print*,'Eos_getD: internalEnergy<eos_smalle!',internalEnergy,eos_smalle
#endif
              internalEnergy = max(internalEnergy, eos_smalle)
              eosData(eint+n) = internalEnergy

              if (eint1_map /= NONEXISTENT) then
                 internalIonEnergy  = solnData(eint1_map,i,j,k)
                 internalIonEnergy = max(internalIonEnergy, eos_smallEion)
              else
                 internalIonEnergy = 0.0
              end if
              eosData(eintIon+n) = internalIonEnergy

              if (eint2_map /= NONEXISTENT) then
                 internalEleEnergy  = solnData(eint2_map,i,j,k)
                 internalEleEnergy = max(internalEleEnergy, eos_smallEele)
              else
                 internalEleEnergy = 0.0
              end if
              eosData(eintEle+n) = internalEleEnergy

              if (eint3_map /= NONEXISTENT) then
                 internalRadEnergy  = solnData(eint3_map,i,j,k)
                 internalRadEnergy = max(internalRadEnergy, eos_smallErad)
              else
                 internalRadEnergy = 0.0
              end if
              eosData(eintRad+n) = internalRadEnergy

              eosData(pres+n) = solnData(pres_map,i,j,k)
              eosData(dens+n) = solnData(dens_map,i,j,k)
              eosData(temp+n) = solnData(temp_map,i,j,k)
              eosData(gamc+n) = solnData(gamc_map,i,j,k)
              if (temp1_map /= NONEXISTENT) then
                 eosData(tempIon+n) = solnData(temp1_map,i,j,k)
              else
                 eosData(tempIon+n) = eosData(temp+n)
              end if
              if (temp2_map /= NONEXISTENT) then
                 eosData(tempEle+n) = solnData(temp2_map,i,j,k)
              else
                 eosData(tempEle+n) = eosData(temp+n)
              end if
              if (temp3_map /= NONEXISTENT) then
                 eosData(tempRad+n) = solnData(temp3_map,i,j,k)
                 if (eosData(tempIon+n)==0.0 .and. eosData(tempEle+n)==0.0 .and. eosData(tempRad+n)==0.0) &
                      eosData(tempRad+n) = eos_smallT
              else
                 eosData(tempRad+n) = eosData(temp+n)
                 if (eosData(tempIon+n)==0.0 .and. eosData(tempEle+n)==0.0 .and. eosData(tempRad+n)==0.0) &
                      eosData(tempRad+n) = eos_smallT
              end if
              if (pres1_map /= NONEXISTENT) then
                 eosData(presIon+n) = solnData(pres1_map,i,j,k)
              else
                 eosData(presIon+n) = eosData(pres+n)
              end if
              if (pres2_map /= NONEXISTENT) then
                 eosData(presEle+n) = solnData(pres2_map,i,j,k)
              else
                 eosData(presEle+n) = eosData(pres+n)
              end if
              if (pres3_map /= NONEXISTENT) then
                 eosData(presRad+n) = solnData(pres3_map,i,j,k)
              else
                 eosData(presRad+n) = eosData(pres+n)
              end if
              if (sele_map /= NONEXISTENT) then
#ifdef EOS_ENTRELE
                 eosData(entrEle+n) = solnData(sele_map,i,j,k)
#else
                 call Driver_abortFlash("EOS_ENTRELE not defined in Eos_getData - using wrong Eos.h?")
#endif
              end if
              if (srad_map /= NONEXISTENT) then
#ifdef EOS_ENTRRAD
                 eosData(entrRad+n) = solnData(srad_map,i,j,k)
#else
                 call Driver_abortFlash("EOS_ENTRRAD not defined in Eos_getData - using wrong Eos.h?")
#endif
              end if
              if(sumy_map /= NONEXISTENT) then
              !! cal says abar=1/sumy
                 if (solnData(sumy_map,i,j,k).NE.0.0) then
                    eosData(abar+n) =  1.0 /  solnData(sumy_map,i,j,k)
                 else
                    eosData(abar+n) =  huge(solnData(sumy_map,i,j,k))
                 end if
              endif
              if((ye_map /= NONEXISTENT).and.(sumy_map /= NONEXISTENT)) then
              !! cal says zbar=ye / sumy and he claims sumy are never zero
                 if (solnData(sumy_map,i,j,k).NE.0.0) then
                    eosData(zbar+n) = solnData(ye_map,i,j,k) /  solnData(sumy_map,i,j,k)
                 else if (solnData(ye_map,i,j,k)==0.0) then
                    eosData(zbar+n) = 0.0
                 else
                    eosData(zbar+n) = huge(solnData(ye_map,i,j,k))
                 end if
              endif
              if(entr_map /= NONEXISTENT) eosData(entr+n) = solnData(entr_map,i,j,k)
           end do
        end do
     end do
  else
     n = 0
     do k = kb,ke
        do j = jb,je
           do i = ib,ie
              kineticEnergy  = 0.5*(solnData(i,j,k,velx_map)**2 + &
                                     solnData(i,j,k,vely_map)**2 + &
                                     solnData(i,j,k,velz_map)**2)
        !! kineticEnergy holds velocity vector information -- 1/2 * Vmag**2
        !! internalEnergy holds eint (directly)  or energyTotal - ekinetic (calculated),
        !!          depending upon eintSwitch
              if(eint_map /= NONEXISTENT) then
                 internalEnergy  = solnData(i,j,k,eint_map)
                 if (solnData(i,j,k,ener_map) > (1.+ eos_eintSwitch)*kineticEnergy) then
                    internalEnergy = solnData(i,j,k,ener_map) - kineticEnergy
                 end if
              else
                 internalEnergy = solnData(i,j,k,ener_map)-kineticEnergy
              endif

              internalEnergy = max(internalEnergy, eos_smalle)

              n=n+1
              eosData(pres+n) = solnData(i,j,k,pres_map)
              eosData(dens+n) = solnData(i,j,k,dens_map)
              eosData(temp+n) = solnData(i,j,k,temp_map)
              eosData(gamc+n) = solnData(i,j,k,gamc_map)
              eosData(eint+n) = internalEnergy
              eosData(ekin+n) = kineticEnergy
              if (temp1_map /= NONEXISTENT) then
                 eosData(tempIon+n) = solnData(i,j,k,temp1_map)
              else
                 eosData(tempIon+n) = eosData(temp+n)
              end if
              if (temp2_map /= NONEXISTENT) then
                 eosData(tempEle+n) = solnData(i,j,k,temp2_map)
              else
                 eosData(tempEle+n) = eosData(temp+n)
              end if
              if (temp3_map /= NONEXISTENT) then
                 eosData(tempRad+n) = solnData(i,j,k,temp3_map)
              else
                 eosData(tempRad+n) = eosData(temp+n)
              end if

              if (eint1_map /= NONEXISTENT) then
                 eosData(eintIon+n) = solnData(i,j,k,eint1_map)
              else
                 eosData(eintIon+n) = eosData(eint+n)
              end if
              if (eint2_map /= NONEXISTENT) then
                 eosData(eintEle+n) = solnData(i,j,k,eint2_map)
              else
                 eosData(eintEle+n) = eosData(eint+n)
              end if
              if (eint3_map /= NONEXISTENT) then
                 eosData(eintRad+n) = solnData(i,j,k,eint3_map)
              else
                 eosData(eintRad+n) = eosData(eint+n)
              end if
              if (pres1_map /= NONEXISTENT) then
                 eosData(presIon+n) = solnData(i,j,k,pres1_map)
              else
                 eosData(presIon+n) = eosData(pres+n)
              end if
              if (pres2_map /= NONEXISTENT) then
                 eosData(presEle+n) = solnData(i,j,k,pres2_map)
              else
                 eosData(presEle+n) = eosData(pres+n)
              end if
              if (pres3_map /= NONEXISTENT) then
                 eosData(presRad+n) = solnData(i,j,k,pres3_map)
              else
                 eosData(presRad+n) = eosData(pres+n)
              end if
              if (sele_map /= NONEXISTENT) then
#ifdef EOS_ENTRELE
                 eosData(entrEle+n) = solnData(i,j,k,sele_map)
#else
                 call Driver_abortFlash("EOS_ENTRELE not defined in Eos_getData - using wrong Eos.h?")
#endif
              end if
              if (srad_map /= NONEXISTENT) then
#ifdef EOS_ENTRRAD
                 eosData(entrRad+n) = solnData(i,j,k,srad_map)
#else
                 call Driver_abortFlash("EOS_ENTRRAD not defined in Eos_getData - using wrong Eos.h?")
#endif
              end if
              if(sumy_map /= NONEXISTENT) then
                 !! cal says abar=1/sumy
                 eosData(abar+n) =  1.0 /  solnData(i,j,k,sumy_map)
              endif
              if((ye_map /= NONEXISTENT).and.(sumy_map /= NONEXISTENT)) then
                 !! cal says zbar=ye / sumy and he claims sumy are never zero
                 eosData(zbar+n) = solnData(i,j,k,ye_map) /  solnData(i,j,k,sumy_map)
              endif
              eosData(entr+n) = solnData(i,j,k,entr_map)
           end do
        end do
     end do

  endif

  return
end subroutine Eos_getData



