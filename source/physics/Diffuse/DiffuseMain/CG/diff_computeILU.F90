!!****if* source/physics/Diffuse/DiffuseMain/CG/diff_computeILU
!!
!!  NAME 
!!
!!  diff_computeILU
!!
!!  SYNOPSIS
!!
!!  call diff_computeILU(AA, JA, IA, UPTR, NZ)
!!
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

!!REORDER(4): solnVec

subroutine diff_computeILU(AA, JA, IA, UPTR, NZ, N)    

  use Timers_interface, ONLY : Timers_start, Timers_stop
  
  implicit none
  
#include "Flash.h"
#include "constants.h"  
  
  integer, intent (IN)     :: NZ
  integer, intent (IN)     :: N 
  real,    intent (INOUT)  :: AA  (NZ)
  integer,    intent (IN)  :: JA  (NZ)
  integer,    intent (IN)  :: IA  (N+1)
  integer,    intent (IN ) :: UPTR(N)
  
  integer :: i, j, jj, k1, k2
  logical :: diagreached  
  integer :: work (N)  


  call Timers_start("diff_computeILU")

 
  work = 0
  do i = 2, N           
     k1 = IA(i)
     k2 = IA(i+1)-1     
     diagreached = .FALSE.
     
     do j=k1, k2
        work(JA(j)) = j        
     end do
     
     do while (JA(k1) .lt. i .and. .not.(diagreached))                 
        AA(K1) = AA(K1) / AA(uptr(JA(K1)))                  
        
        if (AA(uptr(JA(K1))) == 0.0) then
           write(*,*) "zero pivot"
           pause
        end if
        
        do jj = uptr(JA(K1))+1, ia(JA(K1)+1)-1                       
           if (work(JA(jj)) /= 0) then 
              AA(work(JA(jj))) =  AA(work(JA(jj))) - AA(K1)*AA(jj)                                   
           end  if
        end do
        
        if (JA(k1) .eq.uptr(i)-1) then            
           diagreached = .TRUE.           
        else                   
           k1 = k1 + 1           
        end if
        
     end do
     work = 0     
  end do 
  
  
  call Timers_stop("diff_computeILU")
    
  
end subroutine Diff_computeILU
