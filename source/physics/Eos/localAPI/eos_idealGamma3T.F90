!!****if* source/physics/Eos/localAPI/eos_idealGamma3T
!!
!! NAME
!!
!!  eos_idealGamma3T
!!
!! SYNOPSIS
!!
!!  call eos_idealGamma3T(integer(IN) :: mode,
!!                      integer(IN) :: vecLen,
!!                      real(INOUT) :: eosData(vecLen*EOS_NUM),
!!            optional, integer(IN) :: vecBegin,
!!            optional, integer(IN) :: vecEnd,
!!            optional, integer(IN) :: eosType,
!!            optional, integer(IN) :: material,
!!            optional, real(IN)    :: massFrac(vecLen*NSPECIES),
!!      optional,target,logical(IN) :: mask(EOS_VARS+1:EOS_NUM),
!!            optional, integer(IN) :: componentMask(N_EOS_TEMP)    )
!!
!!  This routine implements the gamma law version of the equation of state.
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
!!  eosType :  Type of EOS, should be specified as a constant like EOS_GAM defined
!!             in Eos.h.
!!             This implementation understands
!!               EOS_GAM
!!             and rejects
!!               all other types.
!!             If this optional argument is not present, EOS_GAM is assumed.
!!
!!  material : Indicates to which material the EOS is to be applied,
!!             in a multi-type multi-material context.
!!             Given as an index into the UNK solution vector,
!!             if valid we should have SPECIES_BEGIN <= material <= SPECIES_END.
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
!!
!!
!!
!! SEE ALSO
!! 
!!  Eos.h    defines the variables used.
!!  Eos_wrapped  sets up the data structure.
!!
!!***

subroutine eos_idealGamma3T(mode, vecLen, eosData, vecBegin,vecEnd, &
     eosType, material, massFrac, mask, componentMask)

!==============================================================================

  implicit none
#include "Eos.h"
#include "Flash.h"
  integer, INTENT(in) :: mode, vecLen
  real,INTENT(inout), dimension(EOS_NUM*vecLen) :: eosData 
  integer,optional,INTENT(in) :: vecBegin,vecEnd
  integer,optional,INTENT(in) :: eosType
  integer,optional,INTENT(in) :: material
  logical, optional, INTENT(in),target,dimension(EOS_VARS+1:EOS_NUM) :: mask
  real, optional, INTENT(in),dimension(NSPECIES*vecLen)    :: massFrac
  integer, optional, dimension(N_EOS_TEMP),INTENT(in)::componentMask

  return
end subroutine eos_idealGamma3T
