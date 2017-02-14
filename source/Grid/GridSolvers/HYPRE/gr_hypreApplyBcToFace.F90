!!****if* source/Grid/GridSolvers/HYPRE/gr_hypreApplyBcToFace
!!
!!  NAME
!!
!!  gr_hypreApplyBcToFace
!!
!!  SYNOPSIS
!!
!!  call gr_hypreApplyBcToFace (integer, intent(IN) :: blkLimitsGC(2,MDIM),
!!                              integer, intent(IN) :: blkLimits(2,MDIM),
!!                              integer, intent(IN) :: part,
!!                              integer, intent(IN) :: var,
!!                              integer, intent(IN) :: iFactorB,
!!                              integer, intent(IN) :: bcType,
!!                              integer, intent(IN) :: direction,
!!                              real,    intent(IN) :: bcValue(2),
!!                              real,    intent(IN) :: dt,
!!                              real,    intent(IN) :: theta,
!!                              real,    intent(IN) :: del,
!!                              integer, intent(IN) :: Lower(MDIM),
!!                              real,    intent(IN) :: scalefactor,
!!                              real,    intent(IN) :: faceArea(:,:,:),
!!                              real,    intent(IN) :: solnVec(:,:,:,:))
!!
!!  DESCRIPTION
!!      This routines modifies Matrix A (as in AX=B) or the right hand side B
!!      to account for the applied physical boundary conditions. The
!!      applied boundary conditions might or might not require changes
!!      to the RHS, particularly the way the diffusion operator has been
!!      implemented, OUTFLOW BC ends up as a no-op.
!!
!! ARGUMENTS
!!
!!   blkLimitsGC  : an array that holds the lower and upper indices of the
!!                  section of block with the guard cells. 
!!   blkLimits    : an array that holds the lower and upper indices of the section
!!                  of block without the guard cells.
!!   part         : HYPRE part to which the block belongs to. Determined based
!!                  on its refinement level.
!!   var          : Variable as numbered by HYPRE
!!   iFactorB     : Variable as numbered by FLASH. Conductivity or Opacitiy.
!!   bcType       : GRID_PDE_BND_NEUMANN, VACUUM, OUTSTREAM (in development),
!!                  GRID_PDE_BND_DIRICHLET, GRID_PDE_BND_GIVENGRAD 
!!   direction    : LOWER or UPPER face (i.e. ILO_FACE, IHI_FACE, JLO_FACE ...)
!!   bcValue      : used when GRID_PDE_BND_DIRICHLET/GRID_PDE_BND_GIVENGRAD is specified as bcType
!!                  GRID_PDE_BND_DIRICHLET
!!                    bcValue(1) -> Value of Temperature/Energy-Density specified on boundary.
!!                    bcValue(2) -> Value of Conductivity/Opacity specified on boundary.
!!                  GRID_PDE_BND_GIVENGRAD
!!                     bcValue(1) -> Value of gradient as specified on boundary.
!!                     bcValue(2) -> not used.
!!   dt           : Global time step.
!!   theta        : Varies scheme (0-> Explicit, 1-> backward euler, 0.5 -> Crank Nicholson
!!   del          : dx/dy/dz.
!!   Lower        : Lower limits of the domain (Global, computed based on
!!                  strides and block CornerID).
!!   scalefactor  : used to scale the RHS properly. Used with GIVENRAD, DIRICHLET only.
!!   faceArea     : Areas of the cells on the face (ILO_FACE, JLO_FACE ....).
!!                  At the same refinement level, cell values differ only for non-Cartesian
!!                  coordinates.
!!   solnVec      : slice of UNK corresponding to a particular block (blockID).
!!
!! SIDE EFFECTS
!!
!!  
!! NOTES
!!  
!!
!!***

!!REORDER(4): solnVec

subroutine gr_hypreApplyBcToFace(blkLimits,blkLimitsGC,part,var,iFactorB,bcType,direction, &
     bcValue, dt, theta, del, Lower, scalefactor, faceArea, solnVec, blockID, &
     facBptrX,facBptrY,facBptrZ)
  
  use gr_hypreData,   ONLY : gr_hypreLower, gr_hypreUpper,gr_speedlt,gr_asol, &
                             gr_hypreMatA, gr_hypreVecB
  use Grid_data,      ONLY : gr_meshMe, gr_geometry, gr_meshComm                            
  use Grid_interface,   ONLY : GRID_PDE_BND_PERIODIC,  &
                               GRID_PDE_BND_NEUMANN,   &
                               GRID_PDE_BND_DIRICHLET, &
                               GRID_PDE_BND_GIVENGRAD, &
                               Grid_getCellCoords
  
  implicit none
#include "Flash.h"  
#include "constants.h"
#include "FortranLangFeatures.fh"
#include "HYPREf.h"
  
  integer, intent(IN) :: blkLimits (2,MDIM) 
  integer, intent(IN) :: blkLimitsGC (2,MDIM)
  integer, intent(IN) :: part
  integer, intent(IN) :: var
  integer, intent(IN) :: iFactorB
  integer, intent(IN) :: bcType
  integer, intent(IN) :: direction
  real,    intent(IN) :: bcValue(2)
  real,    intent(IN) :: dt
  real,    intent(IN) :: theta
  real,    intent(IN) :: del
  integer, intent(IN) :: Lower(MDIM)  
  real,    intent(IN) :: scalefactor
  real,    intent(IN) :: faceArea(blkLimitsGC(LOW,IAXIS):blkLimitsGC(HIGH,IAXIS), &
                                  blkLimitsGC(LOW,JAXIS):blkLimitsGC(HIGH,JAXIS), &
                                  blkLimitsGC(LOW,KAXIS):blkLimitsGC(HIGH,KAXIS))
  real,    intent(IN) :: solnVec(NUNK_VARS, blkLimitsGC(LOW,IAXIS):blkLimitsGC(HIGH,IAXIS), &
                                            blkLimitsGC(LOW,JAXIS):blkLimitsGC(HIGH,JAXIS), &
                                            blkLimitsGC(LOW,KAXIS):blkLimitsGC(HIGH,KAXIS))   
  integer, intent(IN) :: blockID
  real, POINTER_INTENT_IN, OPTIONAL :: facBptrX(:,:,:,:),facBptrY(:,:,:,:),facBptrZ(:,:,:,:)
  

  integer :: ierr
  integer :: datasize(MDIM)
  integer :: nentries 
  integer :: stencil_indices(1)  
  integer :: i, j, k
  real    :: coeff  
  real, POINTER, DIMENSION(:,:,:) :: facBptr1
  logical :: useFcB
  integer, dimension(MDIM) :: Upper
  integer :: index    (MDIM)
  integer :: BoxLow   (MDIM)
  integer :: BoxHigh  (MDIM)
  integer :: axis, pos(MDIM), offs(MDIM)
  integer :: ii
  real, parameter :: beta=2./3.
  
  real, allocatable :: BoxVal(:), rCenter(:)
  real, allocatable :: RHSVal(:)
  
  real :: cond,rn,rbdry

  integer :: ioff
  integer :: joff
  integer :: koff

  
  if (bcType == GRID_PDE_BND_NEUMANN) return
      
  BoxLow(IAXIS) = blkLimits(LOW,IAXIS)
  BoxLow(JAXIS) = blkLimits(LOW,JAXIS)
  BoxLow(KAXIS) = blkLimits(LOW,KAXIS)
  
  BoxHigh(IAXIS) = blkLimits(HIGH,IAXIS)
  BoxHigh(JAXIS) = blkLimits(HIGH,JAXIS)
  BoxHigh(KAXIS) = blkLimits(HIGH,KAXIS)

  nullify(facBptr1)
  select case (direction)
  case (ILO_FACE,IHI_FACE)
     if (present(facBptrX)) then
        if (associated(facBptrX)) call AssoMed(facBptr1,facBptrX,iFactorB) ! facBptr1(:,:,:) => facBptrX(iFactorB,:,:,:)
     end if
  case (JLO_FACE,JHI_FACE)     
     if (present(facBptrY)) then
        if (associated(facBptrY)) call AssoMed(facBptr1,facBptrY,iFactorB) ! facBptr1(:,:,:) => facBptrY(iFactorB,:,:,:)
     end if
  case (KLO_FACE,KHI_FACE)
     if (present(facBptrZ)) then
        if (associated(facBptrZ)) call AssoMed(facBptr1,facBptrZ,iFactorB) ! facBptr1(:,:,:) => facBptrZ(iFactorB,:,:,:)
     end if
  end select   
  useFcB = associated(facBptr1)

  index = Lower
  
  ioff = 0
  joff = 0
  koff = 0
  select case (direction)
  case (ILO_FACE)
     BoxHigh(IAXIS) = blkLimits(LOW,IAXIS)
  case (IHI_FACE)
     index  (NDIM-IAXIS+1)  = Lower(NDIM-IAXIS+1)+ &
          (blkLimits(HIGH,IAXIS)-blkLimits(LOW,IAXIS))
     BoxLow (IAXIS)  = blkLimits(HIGH,IAXIS)    
     ioff = 1
  case (JLO_FACE)
     BoxHigh(JAXIS) = blkLimits(LOW,JAXIS)     
  case (JHI_FACE)     
     index  (NDIM-K2D*JAXIS+1)  = Lower(NDIM-K2D*JAXIS+1)+ &
          (blkLimits(HIGH,JAXIS)-blkLimits(LOW,JAXIS))
     BoxLow (JAXIS) = blkLimits(HIGH,JAXIS)
     joff = 1
  case (KLO_FACE)
     BoxHigh(KAXIS) = blkLimits(LOW,KAXIS)
  case (KHI_FACE)
     index  (NDIM-K3D*KAXIS+1)  = Lower(NDIM-K3D*KAXIS+1)+ &
          (blkLimits(HIGH,KAXIS)-blkLimits(LOW,KAXIS))
     BoxLow(KAXIS)  = blkLimits(HIGH,KAXIS)     
     koff = 1
  end select   
  
  do i=1, NDIM
     Upper(i) =  index(i) +  &
          BoxHigh(NDIM-i+1)-BoxLow(NDIM-i+1)
  end do
  
  datasize = BoxHigh - BoxLow + 1
  
  allocate(BoxVal(product(dataSize(1:NDIM))))
  
  select case (bcType)
     
  case (VACUUM)     
     do k = BoxLow(KAXIS), BoxHigh(KAXIS)
        do j = BoxLow(JAXIS), BoxHigh(JAXIS)
           do i = BoxLow(IAXIS), BoxHigh(IAXIS)
              
              ii = (k - BoxLow(KAXIS)  + 1)                             +  &
                   (j - BoxLow(JAXIS))*dataSize(KAXIS)                  +  &
                   (i - BoxLow(IAXIS))*dataSize(KAXIS)*dataSize(JAXIS)   
              
              if (useFcB) then
                 cond  = facBptr1(i+ioff,j+joff,k+koff)* faceArea(i+ioff,j+joff,k+koff)
                 coeff = 2.0*facBptr1(i+ioff,j+joff,k+koff)/(gr_speedlt*del)
              else
                 cond  = solnVec(iFactorB,i,j,k)* faceArea(i+ioff,j+joff,k+koff)
                 coeff = 2.0*solnVec(iFactorB,i,j,k)/(gr_speedlt*del)
              end if
              
              BoxVal(ii) =  - ((coeff - 0.5)/(coeff + 0.5))*cond*theta*dt/(del) & 
                   +  theta*(cond/del)*dt

           end do
        end do
     end do


  case(OUTSTREAM)

     select case (direction)
     case(ILO_FACE,IHI_FACE)
        axis = IAXIS
     case(JLO_FACE,JHI_FACE)
        axis = JAXIS
     case(KLO_FACE,KHI_FACE)
        axis = KAXIS
     end select
     allocate(rCenter(blkLimits(LOW,axis):blkLimits(HIGH,axis)))
     call Grid_getCellCoords(axis,blockID, CENTER, .false.,rCenter,size(rCenter))

     

    do k = BoxLow(KAXIS), BoxHigh(KAXIS)
        do j = BoxLow(JAXIS), BoxHigh(JAXIS)
           do i = BoxLow(IAXIS), BoxHigh(IAXIS)

              pos = (/i,j,k/)
              rn    = rCenter(pos(axis))
              offs = (/ioff,joff,koff/)
              rbdry = rn + del*(offs(axis) - 0.5)   
              
              ii = (k - BoxLow(KAXIS)  + 1)                             +  &
                   (j - BoxLow(JAXIS))*dataSize(KAXIS)                  +  &
                   (i - BoxLow(IAXIS))*dataSize(KAXIS)*dataSize(JAXIS)   
              
              
!              cond  = solnVec(iFactorB,i,j,k)* faceArea(i+ioff,j+joff,k+koff)
!              coeff = 2.0*solnVec(iFactorB,i,j,k)/(gr_speedlt*del)              
              
              BoxVal(ii) =  theta*gr_speedlt*((rn/rbdry)**2)*faceArea(i+ioff,j+joff,k+koff)*dt

           end do
        end do
     end do
     deallocate(rCenter)
     
  case (GRID_PDE_BND_DIRICHLET)    
     allocate(RHSVal(product(dataSize(1:NDIM))))
     RHSVal = 0.0     
     do k = BoxLow(KAXIS), BoxHigh(KAXIS)
        do j = BoxLow(JAXIS), BoxHigh(JAXIS)
           do i = BoxLow(IAXIS), BoxHigh(IAXIS)              
              ii = (k - BoxLow(KAXIS)  + 1)                             +  &
                   (j - BoxLow(JAXIS))*dataSize(KAXIS)                  +  &
                   (i - BoxLow(IAXIS))*dataSize(KAXIS)*dataSize(JAXIS)  
              
              if (bcValue(2) .GE. 0.0) then
                 cond = bcValue(2)*faceArea(i+ioff,j+joff,k+koff)              
              else if (useFcB) then
                 cond = facBptr1(i+ioff,j+joff,k+koff)*faceArea(i+ioff,j+joff,k+koff)
              else
                 cond = solnVec(iFactorB,i,j,k)*faceArea(i+ioff,j+joff,k+koff)                                                
              end if
              
              BoxVal(ii) = cond*theta*dt/del & !! modify matrix
                   + cond*theta*dt/del                                               
              
              !! Modify RHS              
              RHSVal(ii) = scalefactor*cond*(dt/del)*2.0*bcValue(1)              
           end do
        end do
     end do
     
     if (scalefactor /= 0.0) then             
        call HYPRE_SStructVectorAddToBoxValu(gr_hypreVecB, part,index(1:NDIM), &
             Upper(1:NDIM), var, RHSVal(:), ierr)         
     end if
     
     deallocate(RHSVal)


  case (GRID_PDE_BND_GIVENGRAD)
     
     allocate(RHSVal(product(dataSize(1:NDIM))))
     
     RHSVal = 0.0     
     BoxVal = 0.0
     
     do k = BoxLow(KAXIS), BoxHigh(KAXIS)
        do j = BoxLow(JAXIS), BoxHigh(JAXIS)
           do i = BoxLow(IAXIS), BoxHigh(IAXIS)

              ii = (k - BoxLow(KAXIS)  + 1)                             +  &
                   (j - BoxLow(JAXIS))*dataSize(KAXIS)                  +  &
                   (i - BoxLow(IAXIS))*dataSize(KAXIS)*dataSize(JAXIS)             
              
              if (useFcB) then
                 cond = facBptr1(i+ioff,j+joff,k+koff)*faceArea(i+ioff,j+joff,k+koff)
              else
                 cond = solnVec(iFactorB,i,j,k)*faceArea(i+ioff,j+joff,k+koff)
              end if
              
              !! Modify matrix.
              !! BoxVal(ii) = -theta*cond*dt/del +  theta*cond*dt/del
              !! The terms cancels with the diagnol element in matrix.
              
              !! Modify RHS              
              RHSVal(ii) = - scalefactor*theta*dt*bcValue(1)                    
              
           end do
        end do
     end do

     if (scalefactor /= 0.0) then             
        call HYPRE_SStructVectorAddToBoxValu(gr_hypreVecB, part,index(1:NDIM), &
             Upper(1:NDIM), var, RHSVal(:), ierr)         
     end if
     
     
     deallocate(RHSVal)     
     
  end select
  
  
  
  nentries = 1
  stencil_indices(1) = 0  
  
  call HYPRE_SStructMatrixAddToBoxValu(gr_hypreMatA, part, index(1:NDIM), & 
       Upper(1:NDIM), var, nentries, stencil_indices(1), BoxVal(:), ierr)  
  

  deallocate (BoxVal)
  
  
  return
  
contains
  subroutine AssoMed(pp, mm, varNo)
    real,POINTER_INTENT_OUT :: pp(:,:,:)
    real,POINTER_INTENT_IN  :: mm(:,:,:,:)
    integer,intent(in) :: varNo
    call AssoFin(pp,mm(varNo,:,:,:),lbound(mm,1),lbound(mm,2),lbound(mm,3),lbound(mm,4))
  end subroutine AssoMed

  subroutine AssoFin(pp, dd, lb1,lb2,lb3,lb4)
    real,POINTER_INTENT_OUT :: pp(:,:,:)
    integer, intent(in) :: lb1,lb2,lb3,lb4
    real,   intent(in),target :: dd(lb2:,lb3:,lb4:)
    pp => dd
  end subroutine AssoFin
end subroutine gr_hypreApplyBcToFace

