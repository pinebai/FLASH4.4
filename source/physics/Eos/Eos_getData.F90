!!****f* source/physics/Eos/Eos_getData
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
!!             Not used in the current implementation.
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


  implicit none
  
#include "Eos.h"
#include "constants.h"
#include "Flash.h"
  
  integer, intent(in) :: axis, vecLen, gridDataStruct
  integer, dimension(MDIM), intent(in) :: pos
  real, dimension(:),intent(OUT) :: eosData
  real,dimension(:),optional,intent(OUT) :: massFrac
  logical, optional, INTENT(INOUT),dimension(EOS_VARS+1:) :: eosMask     
  real,pointer,dimension(:,:,:,:) :: solnData
  eosData=0.0
  if(present(massFrac))massFrac=0.0
  
  return
end subroutine Eos_getData



