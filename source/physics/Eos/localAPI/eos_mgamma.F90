!!****if* source/physics/Eos/localAPI/eos_mgamma
!!
!! NAME
!!
!!  eos_mgamma
!!
!! SYNOPSIS
!!  
!!  call eos_mgamma(integer(IN) :: mode,
!!                  integer(IN) :: vecLen,
!!                  real(INOUT) :: eosData(EOS_NUM*vecLen),
!!        optional, integer(IN) :: vecBegin,
!!        optional, integer(IN) :: vecEnd,
!!        optional, real(IN)    :: massFrac(vecLen*NSPECIES),
!!        optional, logical(IN) :: mask(EOS_VARS+1:EOS_NUM)    )
!!                 
!!
!! DESCRIPTION
!!
!!  This routine implements the mgamma law version of the 3T equation of state.
!!
!!  ARGUMENTS
!!
!!  mode :    Selects the mode of operation of the Eos unit.
!!             The valid values are MODE_DENS_EI, MODE_DENS_PRES and  
!!             MODE_DENS_TEMP as decribed above.
!!
!!  vecLen   : number of points (cells) for which the eosData array is sized.
!!             If vecBegin and vecEnd are not present, this is also the
!!             number of points (cells) for which EOS computation is to be done.
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
!!  vecEnd   : Index of last cell in eosData to handle.
!!             Can be used to limit operation to a subrange of cells, untested.
!!             If not present, the default is vecLen.
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
!! NOTES
!!
!!  NSPECIES is defined in Flash.h.
!!
!!  EOS_VARS and EOS_NUM  are defined in Eos.h.
!!  Calling funtions should included Eos.h, in order to get the definitions of
!!  Eos-specific constants to be able to populate the eosData and mask arrays.
!!  
!!  MODE_DENS_TEMP, MODE_DENS_EI, MODE_DENS_PRES, etc. are defined in constants.h.
!!
!!  User code should not call this implementation routine directly, but
!!  should call Eos and make sure that the desired Multigamma implementation
!!  is included in the simulation configuration.
!!  All code calling the Eos interface should include a 
!!    use Eos_interface 
!!  statement, preferable with "ONLY" attribute, e.g.,
!!    use Eos_interface, ONLY:  Eos
!!  All routines calling this routine directly should include a 
!!    use eos_localInterface
!!  statement, preferable with "ONLY" attribute, e.g.,
!!    use eos_localInterface, ONLY:  eos_mgamma
!!
!!  For Gamma and Multigamma routines, the entropy and entropy derivatives 
!!  calculations have not been confirmed to be correct.  Use with caution.
!!  In THIS implementation, we do not even TRY to return combined entropies.
!!  Hopefully that will change at some point.
!!
!! SEE ALSO
!! 
!!  Eos.h    defines the variables used.
!!  Eos_wrapped  sets up the data structure.
!!
!!***

subroutine eos_mgamma(mode, vecLen, eosData, vecBegin,vecEnd, massFrac, mask)

!==============================================================================

  implicit none
# include "Eos.h"
# include "Flash.h"
  integer, INTENT(in) :: mode, vecLen
  real,INTENT(inout), dimension(EOS_NUM*vecLen) :: eosData 
  integer,optional,INTENT(in) :: vecBegin,vecEnd
  real, optional, INTENT(in),dimension(NSPECIES*vecLen)    :: massFrac
  logical, optional, INTENT(in),dimension(EOS_VARS+1:EOS_NUM) :: mask

  return
end subroutine eos_mgamma
