!!****if* source/physics/Diffuse/localAPI/diff_computeFluxes
!!
!!  NAME 
!!
!!  diff_computeFluxes
!!
!!  SYNOPSIS
!!
!!  call diff_computeFluxes (integer(IN)           :: blockCount,
!!                           integer(IN)           :: blockList,
!!                           integer(IN)           :: iVar,
!!                              real(IN)           :: iFactorB)
!!
!!
!!  DESCRIPTION 
!!      This routine computes AX, Jacobian free, as we do not have to store AX.
!!      The elements of A are built from FD of general diffusion/conduction equation      
!!
!!      A * dV/dt = d/dx(B*dV/dx) + d/dy(B*dV/dy) + d/dx(B*dV/dz) + C*V + D
!!
!!
!! ARGUMENTS
!!      
!!      
!!
!! SIDE EFFECTS
!!      
!!  
!! NOTES:
!!  
!!
!!***

subroutine diff_computeFluxes (blockCount,blockList,iVar,iFactorB)   

  
  implicit none 
  
  integer,intent(IN) :: blockCount
  integer,intent(IN) :: blockList(blockCount)
  integer,intent(IN) :: iVar 
  integer,intent(IN) :: iFactorB
  

end subroutine
