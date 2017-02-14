!!****if* source/physics/Diffuse/DiffuseMain/CG/diff_applyPC
!!
!!  NAME 
!!
!!  diff_applyPC
!!
!!  SYNOPSIS
!!
!!  call diff_applyPC(Z, R, NZ, N, LU, JA, IA, UPTR, blkLimits, blkLimitsGC)
!!
!!  DESCRIPTION 
!!
!!
!! ARGUMENTS
!!
!!
!! SIDE EFFECTS
!!
!!  
!! NOTES:
!!  
!!
!!***

subroutine diff_applyPC (Z, R, NZ, N, LU, JA, IA, UPTR, blkLimits, blkLimitsGC) 
  
  use Timers_interface, ONLY : Timers_start, Timers_stop
  
  implicit none

#include "Flash.h"
#include "constants.h"

  integer, dimension(LOW:HIGH,MDIM),intent(IN) :: blkLimits, blkLimitsGC
  
  integer, intent (IN) :: NZ
  integer, intent (IN) :: N   
  real,    intent (IN) :: LU(NZ)
  integer,    intent (IN) :: JA(NZ)
  integer,    intent (IN) :: IA(N+1)
  integer,    intent (IN) :: UPTR(N)
  
  real, intent(OUT)    :: Z(blkLimitsGC(LOW,IAXIS):blkLimitsGC(HIGH,IAXIS), &
                            blkLimitsGC(LOW,JAXIS):blkLimitsGC(HIGH,JAXIS), &
                            blkLimitsGC(LOW,KAXIS):blkLimitsGC(HIGH,KAXIS))  
  real, intent(IN)     :: R(blkLimitsGC(LOW,IAXIS):blkLimitsGC(HIGH,IAXIS), &
                            blkLimitsGC(LOW,JAXIS):blkLimitsGC(HIGH,JAXIS), &
                            blkLimitsGC(LOW,KAXIS):blkLimitsGC(HIGH,KAXIS))

  
  real, allocatable,dimension(:) :: S
  integer :: i, j, k,k1,k2    
  integer :: pos_ijk
  integer :: datasize(MDIM)
  
  call Timers_start("diff_applyPC")    
  
  !Z = 0.0  

  datasize  (1:MDIM)=blkLimits  (HIGH,1:MDIM)-blkLimits  (LOW,1:MDIM)+1  
  
  allocate (S(N))
  
  do k = blkLimits(LOW, KAXIS), blkLimits(HIGH, KAXIS)      
     do j = blkLimits(LOW, JAXIS), blkLimits(HIGH, JAXIS)
        do i = blkLimits(LOW, IAXIS), blkLimits(HIGH, IAXIS)
           !! i,j,k => position in matrix            
           pos_ijk = (i-blkLimits(LOW,IAXIS)+1) + datasize(IAXIS)*(j-1-blkLimits(LOW,JAXIS)+1) + &
                datasize(IAXIS)*datasize(JAXIS)*(k-1-blkLimits(LOW,KAXIS)+1)           
           S(pos_ijk) = R(i,j,k)                     
        end do
     end do
  end do
  
  !! with LU decomposition, we need to back/forward Solve.
  !! AX = B
  !! LUX = B, UX = Y
  !! LY = B -> Lower diagnol system.
  !! UX = Y -> Upper diagnol system.    
  
  ! LY = B  
  do i = 2, N
     k1 = IA(i)
     k2 = UPTR(i)-1
     do j=k1, k2
        S(i)=S(i) - LU(j)*S(JA(j))
     end do
  end do

  
  

  ! UX = Y 
  S(N) = S(N) / LU(NZ)   
  
  do i=N-1, 1, -1
     k1 = IA(i+1)-1
     k2 = UPTR(i)+1
     do j=k2, k1
        S(i)=(S(i)- LU(j)*S(JA(j)))
     end do
     
     S(i) = S(i) / LU(UPTR(i))     
  end  do
  

  !! Update solution.
  do k = blkLimits(LOW, KAXIS), blkLimits(HIGH, KAXIS)      
     do j = blkLimits(LOW, JAXIS), blkLimits(HIGH, JAXIS)
        do i = blkLimits(LOW, IAXIS), blkLimits(HIGH, IAXIS)
           !! i,j,k => position in matrix            
           pos_ijk = (i-blkLimits(LOW,IAXIS)+1) + datasize(IAXIS)*(j-1-blkLimits(LOW,JAXIS)+1) + &
                datasize(IAXIS)*datasize(JAXIS)*(k-1-blkLimits(LOW,KAXIS)+1)      
           Z(i,j,k) = S(pos_ijk)           
        end do
     end do
  end do
  
  deallocate (S) 
  
  call Timers_stop("diff_applyPC")
  
end subroutine diff_applyPC
