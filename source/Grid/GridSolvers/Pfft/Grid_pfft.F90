!!****if* source/Grid/GridSolvers/Pfft/Grid_pfft
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
#include "Pfft.h"
#include "constants.h"

  use gr_pfftData, ONLY : pfft_globalLen,pfft_inLen, pfft_outLen, &
                        pfft_midLen,pfft_t1Len, pfft_t2Len, pfft_transformType, &
                        pfft_pclBaseDatType, &
                        pfft_comm, pfft_me, pfft_procGrid,&
                        pfft_work1,pfft_work2,&
                        pfft_trigIaxis,pfft_trigJaxis,&
                        pfft_trigKaxis,pfft_scale, pfft_workSize,pfft_ndim, pfft_usableProc
  use gr_pfftInterface, ONLY : gr_pfftDcftForward, gr_pfftDcftInverse, gr_pfftTranspose
  use Timers_interface, ONLY : Timers_start, Timers_stop

  implicit none
  integer, intent(IN) :: direction
  real, dimension(:),intent(IN) :: inArray
  real, dimension(:), intent(OUT) :: outArray

  integer :: numVec 
  real :: invScale

  if (.not.pfft_usableProc) return
  call Timers_start("PFFT")

  invScale = 1.0


  if(pfft_ndim==1) then  !! treat the one dimensional case
     numVec=1
     if(direction==PFFT_FORWARD) then
        call gr_pfftDcftForward(inArray,outArray,pfft_trigIaxis,&
                              pfft_globalLen(IAXIS),pfft_inLen(IAXIS),&
                              numVec,pfft_transformType(IAXIS),pfft_scale(IAXIS))
     else
        call gr_pfftDcftInverse(inArray,outArray,pfft_trigIaxis,&
                              pfft_globalLen(IAXIS),pfft_inLen(IAXIS),&
                              numVec,pfft_transformType(IAXIS),invScale)
     end if !! end of one dimensional case

  elseif(pfft_ndim==2) then !! begin 2D 

     if(direction==PFFT_FORWARD) then
        call gr_pfftDcftForward(inArray,pfft_work1,pfft_trigIaxis,&
             pfft_globalLen(IAXIS),pfft_inLen(IAXIS),&
             pfft_inLen(JAXIS),pfft_transformType(IAXIS),pfft_scale(IAXIS))

        call gr_pfftTranspose(direction,pfft_pclBaseDatType(IAXIS),pfft_work1,&
             pfft_work2,pfft_t1Len,pfft_midLen,&
             pfft_procGrid(JAXIS),pfft_comm(JAXIS))

        call gr_pfftDcftForward(pfft_work2,outArray,pfft_trigJaxis,&
             pfft_globalLen(JAXIS),pfft_midLen(IAXIS),&
             pfft_midLen(JAXIS),pfft_transformType(JAXIS), pfft_scale(JAXIS))

     else
        call gr_pfftDcftInverse(inArray,pfft_work1,pfft_trigJaxis,&
             pfft_globalLen(JAXIS),pfft_midLen(IAXIS),&
             pfft_midLen(JAXIS),pfft_transformType(JAXIS),invScale)
        call gr_pfftTranspose(direction,pfft_pclBaseDatType(IAXIS),pfft_work1,&
                            pfft_work2,pfft_midLen,pfft_t1Len,&
                            pfft_procGrid(JAXIS),pfft_comm(JAXIS))

        call gr_pfftDcftInverse(pfft_work2,outArray,pfft_trigIaxis,&
             pfft_globalLen(IAXIS),pfft_inLen(IAXIS),&
             pfft_inLen(JAXIS),pfft_transformType(IAXIS),invScale)
     end if  !! end 2D

  else !! begin 3D

     if(direction==PFFT_FORWARD) then  !! begin forward

        numVec = pfft_inLen(JAXIS)*pfft_inLen(KAXIS)

        call gr_pfftDcftForward(inArray,pfft_work1,pfft_trigIaxis,&
                              pfft_globalLen(IAXIS),pfft_inLen(IAXIS),&
                              numVec,pfft_transformType(IAXIS),pfft_scale(IAXIS))
        call gr_pfftTranspose(direction,pfft_pclBaseDatType(IAXIS),pfft_work1,&
                            pfft_work2,pfft_t1Len,pfft_midLen,&
                            pfft_procGrid(JAXIS),pfft_comm(JAXIS))

        numVec=pfft_midLen(JAXIS)*pfft_midLen(KAXIS)
        call gr_pfftDcftForward(pfft_work2,pfft_work1,pfft_trigJaxis,&
                              pfft_globalLen(JAXIS),pfft_midLen(IAXIS),&
                              numVec,pfft_transformType(JAXIS),pfft_scale(JAXIS))
        call gr_pfftTranspose(direction,pfft_pclBaseDatType(JAXIS),pfft_work1,&
                            pfft_work2,pfft_t2Len,pfft_outLen,&
                            pfft_procGrid(KAXIS),pfft_comm(KAXIS))

        numVec=pfft_outLen(JAXIS)*pfft_outLen(KAXIS)
        call gr_pfftDcftForward(pfft_work2,outArray,pfft_trigKaxis,&
                              pfft_globalLen(KAXIS),pfft_outLen(IAXIS),&
                              numVec,pfft_transformType(KAXIS),pfft_scale(KAXIS))

     else !! end forward, begin inverse

        numVec=pfft_outLen(JAXIS)*pfft_outLen(KAXIS)
        call gr_pfftDcftInverse(inArray,pfft_work1,pfft_trigKaxis,&
                              pfft_globalLen(KAXIS),pfft_outLen(IAXIS),&
                              numVec,pfft_transformType(KAXIS),invScale)

        call gr_pfftTranspose(direction,pfft_pclBaseDatType(JAXIS),pfft_work1,&
                            pfft_work2,pfft_outLen,pfft_t2Len,&
                            pfft_procGrid(KAXIS),pfft_comm(KAXIS))

        numVec=pfft_midLen(JAXIS)*pfft_midLen(KAXIS)
        call gr_pfftDcftInverse(pfft_work2,pfft_work1,pfft_trigJaxis,&
                              pfft_globalLen(JAXIS),pfft_midLen(IAXIS),&
                              numVec,pfft_transformType(JAXIS),invScale)
        call gr_pfftTranspose(direction,pfft_pclBaseDatType(IAXIS),pfft_work1,&
                            pfft_work2,pfft_midLen,pfft_t1Len,&
                            pfft_procGrid(JAXIS),pfft_comm(JAXIS))
        numVec=pfft_inLen(JAXIS)*pfft_inLen(KAXIS)
        call gr_pfftDcftInverse(pfft_work2,outArray,pfft_trigIaxis,&
                              pfft_globalLen(IAXIS),pfft_inLen(IAXIS),&
                              numVec,pfft_transformType(IAXIS),invScale)
     end if !! end inverse 3D
  end if !! end 3D

  call Timers_stop("PFFT")
  return
end subroutine Grid_pfft
