!!****if* source/physics/Eos/localAPI/eos_tabIonmix
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
!!  eosType  : the type of eos to be applied
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

subroutine eos_tabIonmix(mode, vecLen, eosData, vecBegin, vecEnd, eosType, subtype, material, mask)

!==============================================================================

  implicit none
# include "Eos.h"
  integer, INTENT(in) :: mode, vecLen
  real,INTENT(inout), dimension(EOS_NUM*vecLen) :: eosData 
  integer,optional,INTENT(in) :: vecBegin,vecEnd
  integer,optional,INTENT(in) :: eosType,subtype
  integer,optional,INTENT(in) :: material
  logical, optional, INTENT(in),target,dimension(EOS_VARS+1:EOS_NUM) :: mask

  return
end subroutine eos_tabIonmix
