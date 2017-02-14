!!****if* source/Grid/GridSolvers/Pfft/gr_pfftIntersperse
!!
!! NAME
!!
!!  gr_pfftIntersperse
!!
!! SYNOPSIS
!!  
!!  gr_pfftIntersperse(real(IN)   :: inArray(:),
!!                     real(OUT)  :: outArray(:),
!!                     integer(IN):: nx,
!!                     integer(IN):: ny,
!!                     integer(IN):: nx1)
!! 
!! DESCRIPTION
!!
!!   This is a helper routine for the forward distributed transpose.
!!   It only moves the data around, but doesn't change it.
!!   If we were to visualize inArray as having the dimensions
!!   (nx/nx1,ny,nx1), then the outArray is shaped as (nx,ny)
!!   
!! ARGUMENTS
!!  
!!  inArray - Array containing input data
!!  outArray - Array for storing transformed data
!!  nx - length of the vector to be split
!!  ny - the undisturbed length
!!  nx1 - the factor by which nx is split
!! 
!! EXAMPLE
!!
!!  For nx=8, ny=4 and nx1=2, the array changes shape as follows
!!
!!   1,1 1,2 1,3 1,4 5,1 5,2 5,3 5,4 
!!   2,1 2,2 2,3 2,4 6,1 6,2 6,3 6,4 
!!   3,1 3,2 3,3 3,4 7,1 7,2 7,3 7,4 
!!   4,1 4,2 4,3 4,4 8,1 8,2 8,3 8,4 
!!   
!!
!!  to 
!!  
!!   1,1 1,2 1,3 1,4 
!!   2,1 2,2 2,3 2,4 
!!   3,1 3,2 3,3 3,4 
!!   4,1 4,2 4,3 4,4 
!!   5,1 5,2 5,3 5,4 
!!   6,1 6,2 6,3 6,4 
!!   7,1 7,2 7,3 7,4 
!!   8,1 8,2 8,3 8,4 
!!   
!!   
!! NOTES  
!!  
!!  Also see gr_pfftDisperse, the action of this routine is
!!  exact opposite of that of gr_pfftDisperse
!!
!!
!!***

subroutine gr_pfftIntersperse(inArray,outArray,nx,ny,nx1)
  implicit none
  integer,intent(IN) :: nx,ny,nx1
  real,dimension(:), intent(IN) :: inArray
  real,dimension(:), intent(OUT) :: outArray
  integer :: j,n,bbeg,bend,abeg,aend,bstart
  logical :: done

  bbeg = 1
  bstart=0
  do j=0,ny-1
     done = .false.
     n=0
     do while(.not.done)
        abeg = 1+nx*(j+ny*n)
        aend = abeg+nx-1
        bend = bbeg+nx-1
        if(bend-bstart >= nx1) then
           aend = abeg+nx1-n*nx-1
           bend = bbeg+nx1-n*nx-1
           done=.true.
        end if
        outArray(bbeg:bend)=inArray(abeg:aend)
        bbeg=bend+1
        n=n+1
     end do
     bstart=bstart+nx1
  end do

  return
end subroutine gr_pfftIntersperse
