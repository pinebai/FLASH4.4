#include "constants.h"
module gr_hypreMultiData
  implicit none

  integer,save,dimension(VARDESC_SIZE) :: gr_dfsvBaseVarDesc

  integer,save :: gr_dfsvPhase
  integer,save :: gr_dfsvFirstUnkVar, gr_dfsvLastUnkVar
  real   ,save :: gr_dfsvTheta,gr_dfsvThetaC,gr_dfsvThetaD
  real   ,save :: gr_dfsvDt

end module gr_hypreMultiData
