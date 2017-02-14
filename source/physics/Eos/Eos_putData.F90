!!****f* source/physics/Eos/Eos_putData
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


subroutine Eos_putData(axis,pos,vecLen,solnData,gridDataStruct,eosData)
  
  implicit none

#include "constants.h"

  integer, intent(in) :: axis,vecLen,gridDataStruct
  integer,dimension(MDIM), intent(in) :: pos
  real,intent(IN) :: eosData(:)
  real, pointer, dimension(:,:,:,:) :: solnData
  
  return
end subroutine Eos_putData



