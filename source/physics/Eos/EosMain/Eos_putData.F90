!!****if* source/physics/Eos/EosMain/Eos_putData
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
  use Eos_data, ONLY : eos_mapLookup, eos_meshMe

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
  integer :: ib,ie,jb,je,kb,ke

  integer :: pres_map,entr_map,gamc_map,temp_map
  integer :: eint_map,ener_map,game_map

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
        print*,'solnData(EINT_VAR,ib:ie,jb,kb)=',solnData(EINT_VAR,ib:ie,jb,kb)
        print*,'eosData (eint+1:eint+vecLen)  =',eosData(eint+1:eint+vecLen)
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

  if(gridDataStruct == SCRATCH) then
     call Driver_abortFlash("Eos_getData : the use of SCRATCH is deprecated")
  end if
  
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
           
           solnData(game_map,i,j,k) = eosData(pres+n)/&
                (eosData(eint+n) *eosData(dens+n)) +1
           
        end do
     end do
  end do
  
  return
end subroutine Eos_putData



