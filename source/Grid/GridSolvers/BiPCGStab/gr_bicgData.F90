!*******************************************************************************

!  Module:      gr_bicgData()

!  Description: Global type and variable declarations for the BiPCGStab solver.

module gr_bicgData

implicit none

#include "constants.h"  ! for MDIM

!===============================================================================


! Ranges of interior indices for blocks.

integer,save :: ili, iui, jli, jui, kli, kui

! Ranges of exterior indices for blocks.

integer,save :: ile, iue, jle, jue, kle, kue

! Boundary conditions, flag to substract mean on the source.
integer,save, dimension(2*MDIM) :: gr_bicgBndTypes
integer,save :: bicg_bnd_cond

! Interpolation for work
integer, parameter :: interp_work_bipcgstab = 2
integer,save, allocatable, dimension(:) :: interp_mask_work_save

! Flag for initialization for preconditioner
logical,save :: precond_init = .true.

! Flag that defines if boundary guard-cells are filled before A*x operations
logical,save :: gcellflg = .true. 

!===============================================================================

end module gr_bicgData
