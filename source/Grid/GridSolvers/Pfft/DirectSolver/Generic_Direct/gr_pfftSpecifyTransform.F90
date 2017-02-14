!!****if* source/Grid/GridSolvers/Pfft/DirectSolver/Generic_Direct/gr_pfftSpecifyTransform
!!
!! NAME
!!
!! gr_pfftSpecifyTransform
!!
!! SYNOPSIS
!!  
!! gr_pfftSpecifyTransform(integer(OUT) :: transformType(MDIM))
!! 
!! DESCRIPTION
!!
!! Initialises the transform type for the particular solver.
!!
!! ARGUMENTS
!!  
!! transformType - Array containing the type of transform required 
!!                 along each axis (e.g. PFFT_REAL or PFFT_COMPLEX)
!!
!! NOTES
!!
!! Called by gr_pfftInit.
!!
!!***

subroutine gr_pfftSpecifyTransform(transformType,baseDatType, bcTypes)

  use Driver_interface, ONLY : Driver_abortFlash
  use Grid_interface, ONLY : GRID_PDE_BND_PERIODIC, GRID_PDE_BND_NEUMANN, &
       GRID_PDE_BND_DIRICHLET, GRID_PDE_BND_GIVENVAL
  use Grid_data, ONLY : gr_meshMe, gr_domainBC


#include "constants.h"
#include "Flash.h"
#include "Pfft.h"
  implicit none
  integer, dimension(MDIM), intent(OUT) :: transformType
  integer, dimension(0:MDIM), intent(OUT), OPTIONAL :: baseDatType
  integer, dimension(2*MDIM), intent(IN),  OPTIONAL :: bcTypes
  logical, save :: firstcall = .true.
  integer :: i

!!$  transformType(IAXIS) = PFFT_REAL
!!$  transformType(JAXIS) = PFFT_REAL
!!$  transformType(KAXIS) = PFFT_COSQ

  if (present(bcTypes)) then
#ifdef DEBUG_ALL
     print*,'bcTypes was present to Specif:',bcTypes
#endif
     do i=1,NDIM
        if (bcTypes(2*i-1) .NE. bcTypes(2*i)) then
           print*,'Unsupported boundary types in direction',i, &
                ', both sides are different from each other!'
           call Driver_abortFlash('gr_pfftSpecifyTransform: invalid PDE boundary type')
        else if (bcTypes(2*i-1) .eq. GRID_PDE_BND_PERIODIC) then
           transformType(i) = PFFT_REAL
        elseif ((bcTypes(2*i-1) .eq. GRID_PDE_BND_NEUMANN)) then
           transformType(i) = PFFT_COSQ
        elseif ((bcTypes(2*i-1) .eq. GRID_PDE_BND_DIRICHLET) .or. &
             (bcTypes(2*i-1) .eq. GRID_PDE_BND_GIVENVAL)) then
           transformType(i) = PFFT_SINQ
        endif
     end do
  else                          !bcTypes is not present, use GridMain BCs

#ifdef DEBUG_ALL
     print*,'gr_domainBC:',gr_domainBC
#endif
  !! **** THIS IS SPECIFIC TO DIRECTSOLVER ROUTINE ****
  !! X direction:
  if (gr_domainBC(LOW,IAXIS) .eq. PERIODIC) then
     transformType(IAXIS) = PFFT_REAL
  elseif ((gr_domainBC(LOW,IAXIS) .eq. OUTFLOW)    .or. &
          (gr_domainBC(LOW,IAXIS) .eq. NOSLIP_INS) .or. &
          (gr_domainBC(LOW,IAXIS) .eq. SLIP_INS)   .or. &
          (gr_domainBC(LOW,IAXIS) .eq. NEUMANN_INS).or. &
          (gr_domainBC(LOW,IAXIS) .eq. INFLOW_INS)) then
     transformType(IAXIS) = PFFT_COSQ
  elseif (gr_domainBC(LOW,IAXIS) .eq. DIRICHLET) then
     transformType(IAXIS) = PFFT_SINQ
  endif

  !! Y direction:
  if (gr_domainBC(LOW,JAXIS) .eq. PERIODIC) then
     transformType(JAXIS) = PFFT_REAL
  elseif ((gr_domainBC(LOW,JAXIS) .eq. OUTFLOW)    .or. &
          (gr_domainBC(LOW,JAXIS) .eq. NOSLIP_INS) .or. &
          (gr_domainBC(LOW,JAXIS) .eq. SLIP_INS)   .or. &
          (gr_domainBC(LOW,JAXIS) .eq. NEUMANN_INS).or. &
          (gr_domainBC(LOW,JAXIS) .eq. INFLOW_INS)) then
     transformType(JAXIS) = PFFT_COSQ
  elseif (gr_domainBC(LOW,JAXIS) .eq. DIRICHLET) then
     transformType(JAXIS) = PFFT_SINQ
  endif

#if NDIM == 3
  !write(*,*) 'in NDIM=3' 
  !! Z direction:
  if (gr_domainBC(LOW,KAXIS) .eq. PERIODIC) then
     transformType(KAXIS) = PFFT_REAL
  elseif ((gr_domainBC(LOW,KAXIS) .eq. OUTFLOW)    .or. &
          (gr_domainBC(LOW,KAXIS) .eq. NOSLIP_INS) .or. &
          (gr_domainBC(LOW,KAXIS) .eq. SLIP_INS)   .or. &
          (gr_domainBC(LOW,KAXIS) .eq. NEUMANN_INS).or. &
          (gr_domainBC(LOW,KAXIS) .eq. INFLOW_INS)) then
     transformType(KAXIS) = PFFT_COSQ
  elseif (gr_domainBC(LOW,KAXIS) .eq. DIRICHLET) then
     transformType(KAXIS) = PFFT_SINQ
  endif
#else
  transformType(KAXIS) = PFFT_TRANSFORM_NONE
#endif

  end if

  if (firstcall) then

  if (gr_meshMe .eq. MASTER_PE) then
     write(*,*) 'Transforms Used on Pfft Solver:'
     select case(transformType(IAXIS))
     case(PFFT_REAL)
        write(*,*) 'X direction: PFFT_REAL'
     case(PFFT_COSQ)
        write(*,*) 'X direction: PFFT_COSQ'
     case(PFFT_SINQ)
        write(*,*) 'X direction: PFFT_SINQ'
     case default
        write(*,*) 'X direction:', transformType(IAXIS)
     end select

     select case(transformType(JAXIS))
     case(PFFT_REAL)
        write(*,*) 'Y direction: PFFT_REAL'
     case(PFFT_COSQ)
        write(*,*) 'Y direction: PFFT_COSQ'
     case(PFFT_SINQ)
        write(*,*) 'Y direction: PFFT_SINQ'
     end select
#if NDIM == 3
     select case(transformType(KAXIS))
     case(PFFT_REAL)
        write(*,*) 'Z direction: PFFT_REAL'
     case(PFFT_COSQ)
        write(*,*) 'Z direction: PFFT_COSQ'
     case(PFFT_SINQ)
        write(*,*) 'Z direction: PFFT_SINQ'
     case default
        write(*,*) 'Z direction:', transformType(KAXIS)
     end select
#endif
  endif

  firstcall = .false.

  endif

  if (present(baseDatType)) then
     baseDatType(0:MDIM) = PFFT_PCLDATA_REAL
  end if

end subroutine gr_pfftSpecifyTransform
