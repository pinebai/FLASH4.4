!!****if* source/physics/Eos/EosMain/multiTemp/Eos_putData
!! NAME
!!
!!  Eos_putData
!! 
!! SYNOPSIS
!!
!!  call Eos_putData(  integer(IN) :: axis,
!!                     integer(IN) :: pos(MDIM),
!!                     integer(IN) :: vecLen,
!!                  real, pointer  :: solnData(:,:,:,:),
!!                     integer(IN) :: gridDataStruct,
!!                     real(IN)    :: eosData(:))
!!
!!
!! DESCRIPTION
!!
!! Eos_putData puts data from an eosData array into a Grid data structure, usually
!! after data in the eosData array have been updated by an Eos call.
!!
!! The Eos_wrapped function is provided for the user's convenience and acts as a simple
!! wrapper to the Eos interface. The Eos interface uses a single, flexible data
!! structure "eosData" to pass the thermodynamic quantities in and out of the
!! function (see Eos). The wrapper hides formation and use of eosData
!! from the users. The wrapper function uses the Eos_putData function to update
!! certain state variables in the relevant section of the block's storage, a vector 
!! at a time. The function can also be used independently to update a vector in a grid block
!! from the values returned by the call to Eos. The arguments axis, pos and vecLen together 
!! specify the relevant vector.
!!
!! If you want to return the derived quantities defined from EOS_VAR+1:EOS_NUM
!! in Eos.h, then you must use the direct interface Eos().
!!
!!  ARGUMENTS 
!!
!!   
!!   axis : the dimension of the vector in the block's storage
!!   pos  : the starting indices of the vector in the block. Note that the
!!          vector has to provide the starting indices for all dimensions
!!   vecLen : the length of the vector
!!   solnData : the solution data for the current block;
!!              various components (variables) of solnData will have been updated
!!              when Eos_putData returns.
!!   gridDataStruct : the relevant grid data structure, on whose data Eos was applied.
!!                    One of CENTER, FACEVAR{X,Y,Z}, GRIDVAR, defined in constants.h .
!!   eosData : the data structure native to Eos unit, in which the computed values 
!!             of the state variables are returned by Eos
!!
!!
!!  EXAMPLE 
!!      if axis = IAXIS, pos(IAXIS)=1,pos(JAXIS)=1,pos(KAXIS)=1 and vecLen=4
!!      then data from applying Eos() to four cells in the first row along IAXIS
!!      of the lower left hand corner of the guard cells in the block is put
!!      into corresponding parts of the Grid data structure.
!!
!!      However if the value were
!!         pos(IAXIS)=iguard+1,
!!         pos(JAXIS)=jguard+1,
!!         pos(KAXIS)=kguard+1, vecLen = NYB, and axis = JAXIS
!!      then data from applying Eos() to the first column along Y axis in the
!!      interior of the block is returned.
!!
!!  NOTES
!!
!!      This interface is called from Eos_wrappped, and is normally not called
!!      by user code directly.
!!
!!      The actual arguments in a call should match those used in a preceding
!!      Eos_getData call used to set up the eosData array.
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
!!     Eos_getData
!!     Eos
!!     Eos.h
!!
!!***

! solnData depends on the ordering on unk
!!REORDER(4): solnData

#define DEBUG_EOS

subroutine Eos_putData(axis,pos,vecLen,solnData,gridDataStruct,eosData)

  use Driver_interface, ONLY : Driver_abortFlash
  use Logfile_interface, ONLY: Logfile_stampMessage 
  use Eos_data, ONLY : eos_meshMe, eos_mapLookup

  implicit none

#include "Eos.h"
#include "constants.h"
#include "Flash.h"
#include "Eos_map.h"

  integer, intent(in) :: axis,vecLen, gridDataStruct
  integer,dimension(MDIM), intent(in) :: pos
  real,intent(IN) :: eosData(:)
  real, pointer:: solnData(:,:,:,:)


  integer :: i,j,k,n, pres,dens,gamc,temp,abar,zbar,eint,entr,ekin
  integer :: tempIon,tempEle,tempRad
  integer :: eintIon,eintEle,eintRad
  integer :: presIon,presEle,presRad
  integer :: ib,ie,jb,je,kb,ke

  integer :: pres_map,entr_map,gamc_map,temp_map
  integer :: eint_map,ener_map,game_map
  integer :: temp1_map,temp2_map,temp3_map
  integer :: eint1_map,eint2_map,eint3_map
  integer :: pres1_map,pres2_map,pres3_map
  integer :: sumy_map,ye_map
  integer :: entrEle, sele_map
  integer :: entrRad, srad_map

  real :: Ye

  integer,allocatable,dimension(:) :: iFlag


  ! check for zero values before calculating gamma
  ! These integers are indexes into the location in eosData just before the storage area for the appropriate variable.
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
  abar = (EOS_ABAR-1)*vecLen
  zbar = (EOS_ZBAR-1)*vecLen
  entr = (EOS_ENTR-1)*vecLen
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
#ifdef EOS_ENTRELE
  entrEle = (EOS_ENTRELE-1)*vecLen
#endif
#ifdef EOS_ENTRRAD
  entrRad = (EOS_ENTRRAD-1)*vecLen
#endif
#ifdef DEBUG_EOS
  allocate(iFlag(vecLen))
  iFlag = 0
  where ( (eosData(eint+1:eint+vecLen) .eq. 0.) .or. (eosData(dens+1:dens+vecLen) .eq. 0.))
     iFlag(1:vecLen) = 1
  end where
  
  !maybe there was a wrong flag set
  if (maxval(iFlag) .gt. 0) then
     if (eos_meshMe .EQ. MASTER_PE) then
        write(*,*) "ERROR After calling Eos, eosData(EOS_EINT) or eosData(EOS_DENS) are zero"
        print*,'iflag=',iflag
#ifdef EINT_VAR
        print*,'solnData(EINT_VAR,ib:ie,jb,kb)=',solnData(EINT_VAR,ib:ie,jb,kb)
#endif
        print*,'eosData (eint+1:eint+vecLen)  =',eosData(eint+1:eint+vecLen)
#ifdef EION_VAR
        print*,'solnData(EION_VAR,ib:ie,jb,kb)=',solnData(EION_VAR,ib:ie,jb,kb)
#endif
        print*,'eosData (eintIon+1:eintIon+vecLen)  =',eosData(eintIon+1:eintIon+vecLen)
#ifdef EELE_VAR
        print*,'solnData(EELE_VAR,ib:ie,jb,kb)=',solnData(EELE_VAR,ib:ie,jb,kb)
#endif
        print*,'eosData (eintEle+1:eintEle+vecLen)  =',eosData(eintEle+1:eintEle+vecLen)
#ifdef ERAD_VAR
        print*,'solnData(ERAD_VAR,ib:ie,jb,kb)=',solnData(ERAD_VAR,ib:ie,jb,kb)
#endif
        print*,'eosData (eintRad+1:eintRad+vecLen)  =',eosData(eintRad+1:eintRad+vecLen)
        print*,'solnData(DENS_VAR,ib:ie,jb,kb)=',solnData(DENS_VAR,ib:ie,jb,kb)
        print*,'eosData(dens+1:dens+vecLen)   =',eosData(dens+1:dens+vecLen)
        write(*,*) "  Perhaps the initialization routine is wrong..... or"
        write(*,*) "  perhaps the runtime parameter eosMode is wrong."
        write(*,*) "     Check constants.h to determine value of MODE_DENS_??"
     endif
     call Logfile_stampMessage('[Eos_putData] ERROR Density or Internal Energy are zero after a call to EOS!')
     call Driver_abortFlash('[Eos_putData] ERROR Density or Internal Energy are zero after a call to EOS!')
  end if
  deallocate(iFlag)
#endif


  ! Initializations:   grab the solution data from UNK and determine
  !   the length of the data being operated upon
  pres_map = eos_mapLookup(EOSMAP_PRES,EOS_OUT,gridDataStruct)
  temp_map = eos_mapLookup(EOSMAP_TEMP,EOS_OUT,gridDataStruct)
  gamc_map = eos_mapLookup(EOSMAP_GAMC,EOS_OUT,gridDataStruct)
  game_map = eos_mapLookup(EOSMAP_GAME,EOS_OUT,gridDataStruct)
  eint_map = eos_mapLookup(EOSMAP_EINT,EOS_OUT,gridDataStruct)
  ener_map = eos_mapLookup(EOSMAP_ENER,EOS_OUT,gridDataStruct)
  entr_map = eos_mapLookup(EOSMAP_ENTR,EOS_OUT,gridDataStruct)
  temp1_map = eos_mapLookup(EOSMAP_TEMP1,EOS_OUT,gridDataStruct)
  temp2_map = eos_mapLookup(EOSMAP_TEMP2,EOS_OUT,gridDataStruct)
  temp3_map = eos_mapLookup(EOSMAP_TEMP3,EOS_OUT,gridDataStruct)
  eint1_map = eos_mapLookup(EOSMAP_EINT1,EOS_OUT,gridDataStruct)
  eint2_map = eos_mapLookup(EOSMAP_EINT2,EOS_OUT,gridDataStruct)
  eint3_map = eos_mapLookup(EOSMAP_EINT3,EOS_OUT,gridDataStruct)
  pres1_map = eos_mapLookup(EOSMAP_PRES1,EOS_OUT,gridDataStruct)
  pres2_map = eos_mapLookup(EOSMAP_PRES2,EOS_OUT,gridDataStruct)
  pres3_map = eos_mapLookup(EOSMAP_PRES3,EOS_OUT,gridDataStruct)
  sele_map = eos_mapLookup(EOSMAP_SELE,EOS_OUT,gridDataStruct)
  srad_map = eos_mapLookup(EOSMAP_SRAD,EOS_OUT,gridDataStruct)

  sumy_map = eos_mapLookup(EOSMAP_SUMY,EOS_OUT,gridDataStruct)
  ye_map = eos_mapLookup(EOSMAP_YE,EOS_OUT,gridDataStruct)


  if(gridDataStruct /= SCRATCH) then
     n=0
     do k = kb,ke
        do j = jb,je
           do i = ib,ie
              n=n+1
              solnData(pres_map,i,j,k) = eosData(pres+n)
              solnData(temp_map,i,j,k) = eosData(temp+n)
              solnData(gamc_map,i,j,k) = eosData(gamc+n)
              if(eint_map /= NONEXISTENT)solnData(eint_map,i,j,k) = eosData(eint+n)
              if(ener_map /= NONEXISTENT)solnData(ener_map,i,j,k) = eosData(eint+n)+eosData(ekin+n) 
              if(entr_map /= NONEXISTENT)solnData(entr_map,i,j,k) = eosData(entr+n)
              if(temp1_map /= NONEXISTENT)solnData(temp1_map,i,j,k) = eosData(tempIon+n)
              if(temp2_map /= NONEXISTENT)solnData(temp2_map,i,j,k) = eosData(tempEle+n)
              if(temp3_map /= NONEXISTENT)solnData(temp3_map,i,j,k) = eosData(tempRad+n)
              if(eint1_map /= NONEXISTENT)solnData(eint1_map,i,j,k) = eosData(eintIon+n)
              if(eint2_map /= NONEXISTENT)solnData(eint2_map,i,j,k) = eosData(eintEle+n)
              if(eint3_map /= NONEXISTENT)solnData(eint3_map,i,j,k) = eosData(eintRad+n)
              if(pres1_map /= NONEXISTENT)solnData(pres1_map,i,j,k) = eosData(presIon+n)
              if(pres2_map /= NONEXISTENT)solnData(pres2_map,i,j,k) = eosData(presEle+n)
              if(pres3_map /= NONEXISTENT)solnData(pres3_map,i,j,k) = eosData(presRad+n)
              if(ye_map /= NONEXISTENT) then
                 Ye = eosData(zbar+n) / eosData(abar+n)
                 solnData(ye_map,i,j,k) = Ye
#ifdef ZBAR_VAR
                 solnData(ZBAR_VAR,i,j,k) = eosData(zbar+n)
#endif
              end if
              if(sumy_map /= NONEXISTENT) then
                 solnData(sumy_map,i,j,k) = 1.0/eosData(abar+n)
              end if
              if(sele_map /= NONEXISTENT)then
#ifdef EOS_ENTRELE
                 solnData(sele_map,i,j,k) = eosData(entrEle+n)
#else
                 call Driver_abortFlash("EOS_ENTRELE not defined in Eos_putData - using wrong Eos.h?")
#endif
              end if
              if(srad_map /= NONEXISTENT)then
#ifdef EOS_ENTRRAD
                 solnData(srad_map,i,j,k) = eosData(entrRad+n)
#else
                 call Driver_abortFlash("EOS_ENTRRAD not defined in Eos_putData - using wrong Eos.h?")
#endif
              end if

              solnData(game_map,i,j,k) = eosData(pres+n)/&
                   (eosData(eint+n) *eosData(dens+n)) +1


#if defined(TITE_VAR) && defined(TION_VAR) && defined(TELE_VAR)
              if (solnData(TELE_VAR,i,j,k) .NE. 0.0) then
                 solnData(TITE_VAR,i,j,k) = solnData(TION_VAR,i,j,k) / solnData(TELE_VAR,i,j,k)
              else if (solnData(TION_VAR,i,j,k) == 0.0) then
                 solnData(TITE_VAR,i,j,k) = -1.0
              else
                 solnData(TITE_VAR,i,j,k) = sign(huge(solnData(TITE_VAR,i,j,k)),solnData(TION_VAR,i,j,k))
              end if
#endif
#if defined(PIPE_VAR) && defined(PION_VAR) && defined(PELE_VAR)
              if (solnData(PELE_VAR,i,j,k) .NE. 0.0) then
                 solnData(PIPE_VAR,i,j,k) = solnData(PION_VAR,i,j,k) / solnData(PELE_VAR,i,j,k)
              else if (solnData(PION_VAR,i,j,k) == 0.0) then
                 solnData(PIPE_VAR,i,j,k) = -6.0
              else
                 solnData(PIPE_VAR,i,j,k) = sign(huge(solnData(PIPE_VAR,i,j,k)),solnData(PION_VAR,i,j,k))
              end if
#endif

#if defined(TRTE_VAR) && defined(TRAD_VAR) && defined(TELE_VAR)
              if (solnData(TELE_VAR,i,j,k) .NE. 0.0) then
                 solnData(TRTE_VAR,i,j,k) = solnData(TRAD_VAR,i,j,k) / solnData(TELE_VAR,i,j,k)
              else if (solnData(TRAD_VAR,i,j,k) == 0.0) then
                 solnData(TRTE_VAR,i,j,k) = -1.0
              else
                 solnData(TRTE_VAR,i,j,k) = sign(huge(solnData(TRTE_VAR,i,j,k)),solnData(TRAD_VAR,i,j,k))
              end if
#endif
#if defined(PRPE_VAR) && defined(PRAD_VAR) && defined(PELE_VAR)
              if (solnData(PELE_VAR,i,j,k) .NE. 0.0) then
                 solnData(PRPE_VAR,i,j,k) = solnData(PRAD_VAR,i,j,k) / solnData(PELE_VAR,i,j,k)
              else if (solnData(PRAD_VAR,i,j,k) == 0.0) then
                 solnData(PRPE_VAR,i,j,k) = -6.0
              else
                 solnData(PRPE_VAR,i,j,k) = sign(huge(solnData(PRPE_VAR,i,j,k)),solnData(PRAD_VAR,i,j,k))
              end if
#endif

           end do
        end do
     end do
  else
     n=0
     do k = kb,ke
        do j = jb,je
           do i = ib,ie
              n=n+1
              solnData(i,j,k,pres_map) = eosData(pres+n)
              solnData(i,j,k,temp_map) = eosData(temp+n)
              solnData(i,j,k,gamc_map) = eosData(gamc+n)
              if(eint_map /= NONEXISTENT)solnData(i,j,k,eint_map) = eosData(eint+n)
              if(ener_map /= NONEXISTENT)solnData(i,j,k,ener_map) = eosData(eint+n)+eosData(ekin+n) 
              if(entr_map /= NONEXISTENT)solnData(i,j,k,entr_map) = eosData(entr+n)
              if(temp1_map /= NONEXISTENT)solnData(i,j,k,temp1_map) = eosData(tempIon+n)
              if(temp2_map /= NONEXISTENT)solnData(i,j,k,temp2_map) = eosData(tempEle+n)
              if(temp3_map /= NONEXISTENT)solnData(i,j,k,temp3_map) = eosData(tempRad+n)
              if(eint1_map /= NONEXISTENT)solnData(i,j,k,eint1_map) = eosData(eintIon+n)
              if(eint2_map /= NONEXISTENT)solnData(i,j,k,eint2_map) = eosData(eintEle+n)
              if(eint3_map /= NONEXISTENT)solnData(i,j,k,eint3_map) = eosData(eintRad+n)
              if(pres1_map /= NONEXISTENT)solnData(i,j,k,pres1_map) = eosData(presIon+n)
              if(pres2_map /= NONEXISTENT)solnData(i,j,k,pres2_map) = eosData(presEle+n)
              if(pres3_map /= NONEXISTENT)solnData(i,j,k,pres3_map) = eosData(presRad+n)
              if(ye_map /= NONEXISTENT) then
                 Ye = eosData(zbar+n) / eosData(abar+n)
                 solnData(i,j,k,ye_map) = Ye
              end if
              if(sumy_map /= NONEXISTENT) then
                 solnData(i,j,k,sumy_map) = 1.0/eosData(abar+n)
              end if
              if(sele_map /= NONEXISTENT)then
#ifdef EOS_ENTRELE
                 solnData(i,j,k,sele_map) = eosData(entrEle+n)
#else
                 call Driver_abortFlash("EOS_ENTRELE not defined in Eos_putData - using wrong Eos.h?")
#endif
              end if
              if(srad_map /= NONEXISTENT)then
#ifdef EOS_ENTRRAD
                 solnData(i,j,k,srad_map) = eosData(entrRad+n)
#else
                 call Driver_abortFlash("EOS_ENTRRAD not defined in Eos_putData - using wrong Eos.h?")
#endif
              end if
              solnData(i,j,k,game_map) = eosData(pres+n)/&
                   (eosData(eint+n) *eosData(dens+n)) +1
           end do
        end do
     end do

  end if
  
  return
end subroutine Eos_putData



