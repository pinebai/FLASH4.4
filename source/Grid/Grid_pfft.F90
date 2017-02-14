!!****f* source/Grid/Grid_pfft
!!
!! NAME
!!
!!  Grid_pfft
!!
!! SYNOPSIS
!!
!!  call Grid_pfft (integer(IN) :: direction,
!!                  real(IN)    :: inArray(:),
!!                  real(OUT)   :: outArray(:))
!!
!! DESCRIPTION
!!    
!!   This routine is the main interface for the parallel fft solver.
!!   It can do real-to-complex, complex-complex, sine and cosine 
!!   transforms in parallel. The calling routine must first call the
!!   initialization routine, gr_pfftInit to create the data structures
!!   and trignometric tables for calculating the FFT. Once the calling
!!   routine is done, it should call gr_pfftFinalize to deallocate and
!!   clean up. Note that in the same simulation, it is possible to 
!!   to call this routine for different sizes and different types 
!!   of transforms, but each time the sequence must be 
!!    call Grid_pfftInit(....)
!!    call Grid_pfft(PFFT_FORWARD ...
!!     .
!!     .
!!    call Grid_pfft(PFFT_INVERSE ...
!!    call Grid_pfftFinalize()
!!   The routine assumes that the domain is divided into pencils, that
!!   is the first dimension is all within the processor, while the 
!!   remaining dimensions are distributed over a one or two dimensional
!!   processor grid. Note that the data distribution at the end of the
!!   forward transform is different from the original distriution, but
!!   it is restored when the inverse transform is performed.
!!   For example, if the global data is of the size NX,NY,NZ distributed
!!   over py,pz processors, then the original distribution is
!!   NX,NY/py,NZ/pz on each processor. After forward transform the 
!!   distribution is NZ,NX/py,NY/pz. And the original distribution is
!!   restored upon inverse transform.
!!
!! ARGUMENTS
!!  
!!   direction - to indicate whether it is forward or inverse transform
!!   inArray   - single dimension array containing the data to be
!!               transformed
!!   outArray   - array containing transformed data.
!!
!! NOTES 
!!  
!!   Any routine calling this one must include the following line at the top
!! #include "Pfft.h"
!!   Also, please see Grid_pfftInit.
!!
!!*** 

subroutine Grid_pfft(direction,inArray,outArray)
  implicit none
  integer, intent(IN) :: direction
  real, intent(IN),dimension(:) :: inArray
  real, intent(OUT),dimension(:) :: outArray
  outArray = 0.0
  return
end subroutine Grid_pfft
