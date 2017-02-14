!!****if* source/Grid/GridSolvers/Pfft/HomBcTrigSolver/gr_pfftSpecifyTransform
!!
!! NAME
!!
!! gr_pfftSpecifyTransform
!!
!! SYNOPSIS
!!  
!! call gr_pfftSpecifyTransform(integer(OUT) :: transformType(MDIM),
!!                   OPTIONAL,  integer(OUT) :: baseDatType(0:MDIM),
!!                   OPTIONAL,  integer(IN)  :: bcTypes    (2*MDIM))
!! 
!! DESCRIPTION
!!
!! Initialises the transform type for the particular solver.
!!
!! ARGUMENTS
!!  
!! transformType - Array containing the type of transform required 
!!                 along each axis (e.g. PFFT_REAL2C or PFFT_COMPLEX)
!! baseDatType   - If present, indicates the data type as which the
!!                 data after the various stages of transformation
!!                 should be interpreted; should be PFFT_PCLDATA_REAL
!!                 in all elements for this implementation.
!! bcTypes       - If present, boundary conditions to be applied in
!!                 various directions; default is GRID_PDE_BND_PERIODIC
!!                 for all boundaries. Each element should be one of
!!                   GRID_PDE_BND_PERIODIC,
!!                   GRID_PDE_BND_DIRICHLET,
!!                   GRID_PDE_BND_NEUMANN.
!!   
!! NOTES
!!
!! Called by gr_pfftInit.
!!
!!***

subroutine gr_pfftSpecifyTransform(transformType,baseDatType, bcTypes)
  use Grid_interface, ONLY : GRID_PDE_BND_PERIODIC, GRID_PDE_BND_NEUMANN, &
       GRID_PDE_BND_DIRICHLET
  use Driver_interface, ONLY : Driver_abortFlash
#include "constants.h"
#include "Flash.h"
#include "Pfft.h"
  implicit none
  integer, dimension(MDIM), intent(OUT) :: transformType
  integer, dimension(0:MDIM), intent(OUT), OPTIONAL :: baseDatType
  integer, dimension(2*MDIM), intent(IN),  OPTIONAL :: bcTypes

  integer :: i

  if (present(bcTypes)) then
     if (bcTypes(1) /= bcTypes(2)) then
        if ((bcTypes(1)==GRID_PDE_BND_NEUMANN.OR.bcTypes(1)==GRID_PDE_BND_DIRICHLET) .NEQV. &
            (bcTypes(2)==GRID_PDE_BND_NEUMANN.OR.bcTypes(2)==GRID_PDE_BND_DIRICHLET)) &
            call Driver_abortFlash("This Poisson solver requires the same type of boundaries left and right!")
     end if
     if (bcTypes(1) /= GRID_PDE_BND_PERIODIC .AND. bcTypes(1) /= GRID_PDE_BND_NEUMANN &
          .AND. bcTypes(1) /= GRID_PDE_BND_DIRICHLET) &
          call Driver_abortFlash(&
          "This Poisson solver requires periodic or homogeneous Dirichlet or Neumann boundaries in the X direction!")
#if NDIM > 1
     if (bcTypes(3) /= bcTypes(4)) then
        if ((bcTypes(3)==GRID_PDE_BND_NEUMANN.OR.bcTypes(3)==GRID_PDE_BND_DIRICHLET) .NEQV. &
            (bcTypes(4)==GRID_PDE_BND_NEUMANN.OR.bcTypes(4)==GRID_PDE_BND_DIRICHLET)) &
            call Driver_abortFlash("This Poisson solver requires the same type of boundaries up and down!")
     end if
     if (bcTypes(3) /= GRID_PDE_BND_PERIODIC .AND. bcTypes(3) /= GRID_PDE_BND_NEUMANN &
          .AND. bcTypes(3) /= GRID_PDE_BND_DIRICHLET) &
          call Driver_abortFlash(&
          "This Poisson solver requires periodic or homogeneous Dirichlet or Neumann boundaries in the Y direction!")
#endif
#if NDIM > 2
     if (bcTypes(5) /= bcTypes(6)) then
        if ((bcTypes(5)==GRID_PDE_BND_NEUMANN.OR.bcTypes(5)==GRID_PDE_BND_DIRICHLET) .NEQV. &
            (bcTypes(6)==GRID_PDE_BND_NEUMANN.OR.bcTypes(6)==GRID_PDE_BND_DIRICHLET)) &
            call Driver_abortFlash("This Poisson solver requires the same type of boundaries front and back!")
     end if
     if (bcTypes(5) /= GRID_PDE_BND_PERIODIC .AND. bcTypes(5) /= GRID_PDE_BND_NEUMANN &
          .AND. bcTypes(5) /= GRID_PDE_BND_DIRICHLET) &
          call Driver_abortFlash(&
          "This Poisson solver requires periodic or homogeneous Dirichlet or Neumann boundaries in the Z direction!")
#endif

     do i=1,NDIM
        select case (bcTypes(2*i-1))
        case(GRID_PDE_BND_PERIODIC)
           transformType(i) = PFFT_REAL
        case(GRID_PDE_BND_NEUMANN)
           if (bcTypes(2*i)==GRID_PDE_BND_DIRICHLET) then
              transformType(i) = PFFT_COS_IV
           else
              transformType(i) = PFFT_COS_CC
           end if
        case(GRID_PDE_BND_DIRICHLET)
           if (bcTypes(2*i)==GRID_PDE_BND_NEUMANN) then
           transformType(i) = PFFT_SIN_IV
        else
           transformType(i) = PFFT_SIN_CC
        end if
        case default
           call Driver_abortFlash("This Poisson solver requires periodic or &
                &homogeneous Dirichlet or Neumann boundaries!")
        end select
     end do

  else
     transformType(IAXIS) = PFFT_REAL
     transformType(JAXIS:KAXIS) = PFFT_REAL
  end if

  if (present(baseDatType)) then
     baseDatType(0:MDIM) = PFFT_PCLDATA_REAL
  end if
end subroutine gr_pfftSpecifyTransform
