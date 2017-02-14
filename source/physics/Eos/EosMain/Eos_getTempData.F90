!!****if* source/physics/Eos/EosMain/Eos_getTempData
!! NAME
!!
!!  Eos_getTempData
!! 
!! SYNOPSIS
!!
!!  call Eos_getTempData(  integer(IN) :: axis,
!!                     integer(IN) :: pos(MDIM),
!!                     integer(IN) :: vecLen,
!!                  real, pointer  :: solnData(:,:,:,:),
!!                     integer(IN) :: gridDataStruct,
!!                     real(OUT)   :: eosData(:))
!!
!!
!!
!! DESCRIPTION
!!
!! Eos_getTempData gets temperatue data from a Grid data structure into an eosData array, for
!! passing to a subsequent Eos call.
!!
!!
!! While Eos does not know anything about blocks, Eos_getTempData takes its
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


subroutine Eos_getTempData(axis,pos,vecLen,solnData,gridDataStruct,eosData,mode)

  use Eos_data, ONLY: eos_smallT, eos_mapLookup
  use Driver_interface, ONLY : Driver_abortFlash
  implicit none
  
#include "Eos.h"
#include "Eos_map.h"
#include "constants.h"
#include "Flash.h"
  
  integer, intent(in) :: axis, vecLen, gridDataStruct, mode
  integer, dimension(MDIM), intent(in) :: pos
  real, dimension(:),intent(OUT) :: eosData
  real, pointer:: solnData(:,:,:,:)


  integer :: i,j,k,n,m,temp
  integer :: tempIon,tempEle,tempRad
  integer :: temp_map
  integer :: temp1_map,temp2_map,temp3_map
  integer :: ib,ie,jb,je,kb,ke
!! ---------------------------------------------------------------------------------

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
  temp = (EOS_TEMP-1)*vecLen    ! EOS_TEMP is always defined in Eos.h
#ifdef EOS_TEMPION
  tempIon = (EOS_TEMPION-1)*vecLen
#else
  tempIon = temp
#endif
#ifdef EOS_TEMPELE
  tempEle = (EOS_TEMPELE-1)*vecLen
#else
  tempEle = temp
#endif
#ifdef EOS_TEMPRAD
  tempRad = (EOS_TEMPRAD-1)*vecLen
#else
  tempRad = temp
#endif

  temp_map = eos_mapLookup(EOSMAP_TEMP,EOS_IN,gridDataStruct)
  if(mode .NE. MODE_DENS_TEMP) then
     temp1_map = eos_mapLookup(EOSMAP_TEMP1,EOS_IN,gridDataStruct)
     temp2_map = eos_mapLookup(EOSMAP_TEMP2,EOS_IN,gridDataStruct)
     temp3_map = eos_mapLookup(EOSMAP_TEMP3,EOS_IN,gridDataStruct)
  end if

  if(gridDataStruct == SCRATCH) then
     call Driver_abortFlash("Eos_getTempData : the use of SCRATCH is invalid here")
  end if

  n = 0
  do k = kb,ke
     do j = jb,je
        do i = ib,ie
           n=n+1

           eosData(temp+n) = solnData(temp_map,i,j,k)

           if(mode .NE. MODE_DENS_TEMP) then
              if (temp1_map /= NONEXISTENT .AND. tempIon .NE. temp) then
                 eosData(tempIon+n) = solnData(temp1_map,i,j,k)
              else
                 eosData(tempIon+n) = eosData(temp+n)
              end if
              if (temp2_map /= NONEXISTENT .AND. tempEle .NE. temp) then
                 eosData(tempEle+n) = solnData(temp2_map,i,j,k)
              else
                 eosData(tempEle+n) = eosData(temp+n)
              end if
              if (temp3_map /= NONEXISTENT .AND. tempRad .NE. temp) then
                 eosData(tempRad+n) = solnData(temp3_map,i,j,k)
!!$                 if (eosData(tempIon+n)==0.0 .and. eosData(tempEle+n)==0.0 .and. eosData(tempRad+n)==0.0) &
!!$                      eosData(tempRad+n) = eos_smallT
              else
                 eosData(tempRad+n) = eosData(temp+n)
!!$                 if (eosData(tempIon+n)==0.0 .and. eosData(tempEle+n)==0.0 .and. eosData(tempRad+n)==0.0) &
!!$                      eosData(tempRad+n) = eos_smallT
              end if
           end if
        end do
     end do
  end do
  
  return
end subroutine Eos_getTempData 


subroutine Eos_getTempDataFromVec(solnVec,eosData,mode)

  use Eos_data, ONLY: eos_smallT, eos_mapLookup
  implicit none
  
  integer, intent(in) :: mode
  real, dimension(:),intent(OUT) :: eosData
  real, dimension(NUNK_VARS),intent(IN) :: solnVec

  integer,parameter :: vecLen = 1
  integer,parameter :: gridDataStruct = CENTER
  integer :: temp
  integer :: tempIon,tempEle,tempRad
  integer :: temp_map
  integer :: temp1_map,temp2_map,temp3_map
!! ---------------------------------------------------------------------------------

  ! These integers are indexes into the location in eosData just before the storage area for the appropriate variable.
  temp = (EOS_TEMP-1)*vecLen    ! EOS_TEMP is always defined in Eos.h
#ifdef EOS_TEMPION
  tempIon = (EOS_TEMPION-1)*vecLen
#else
  tempIon = temp
#endif
#ifdef EOS_TEMPELE
  tempEle = (EOS_TEMPELE-1)*vecLen
#else
  tempEle = temp
#endif
#ifdef EOS_TEMPRAD
  tempRad = (EOS_TEMPRAD-1)*vecLen
#else
  tempRad = temp
#endif

  temp_map = eos_mapLookup(EOSMAP_TEMP,EOS_IN,gridDataStruct)
  if(mode .NE. MODE_DENS_TEMP) then
     temp1_map = eos_mapLookup(EOSMAP_TEMP1,EOS_IN,gridDataStruct)
     temp2_map = eos_mapLookup(EOSMAP_TEMP2,EOS_IN,gridDataStruct)
     temp3_map = eos_mapLookup(EOSMAP_TEMP3,EOS_IN,gridDataStruct)
  end if

  
  eosData(temp+1) = solnVec(temp_map)

  if(mode .NE. MODE_DENS_TEMP) then
     if (temp1_map /= NONEXISTENT .AND. tempIon .NE. temp) then
        eosData(tempIon+1) = solnVec(temp1_map)
     else
        eosData(tempIon+1) = eosData(temp+1)
     end if
     if (temp2_map /= NONEXISTENT .AND. tempEle .NE. temp) then
        eosData(tempEle+1) = solnVec(temp2_map)
     else
        eosData(tempEle+1) = eosData(temp+1)
     end if
     if (temp3_map /= NONEXISTENT .AND. tempRad .NE. temp) then
        eosData(tempRad+1) = solnVec(temp3_map)
!!$                 if (eosData(tempIon+1)==0.0 .and. eosData(tempEle+1)==0.0 .and. eosData(tempRad+1)==0.0) &
!!$                      eosData(tempRad+1) = eos_smallT
     else
        eosData(tempRad+1) = eosData(temp+1)
!!$                 if (eosData(tempIon+1)==0.0 .and. eosData(tempEle+1)==0.0 .and. eosData(tempRad+1)==0.0) &
!!$                      eosData(tempRad+1) = eos_smallT
     end if
  end if
  
  return
end subroutine Eos_getTempDataFromVec



subroutine Eos_getTempDataF(axis,pos,vecLen,solnData,gridDataStruct,eosData,mode)
  use Eos_interface, ONLY: Eos_getTempData
  implicit none
  integer, intent(in) :: axis, vecLen, gridDataStruct, mode
  integer, dimension(MDIM), intent(in) :: pos
  real, dimension(EOS_NUM*vecLen),intent(OUT) :: eosData
  real, pointer:: solnData(:,:,:,:)

  call  Eos_getTempData(axis,pos,vecLen,solnData,gridDataStruct,eosData,mode)
end subroutine Eos_getTempDataF
