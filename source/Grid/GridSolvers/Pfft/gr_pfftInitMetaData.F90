!!****if* source/Grid/GridSolvers/Pfft/gr_pfftInitMetaData
!!
!! NAME 
!!
!!   gr_pfftInitMetaData
!!
!! SYNOPSIS
!!
!!   gr_pfftInitMetaData(integer(IN) :: ndim)
!!
!! DESCRIPTION 
!!
!!  allocate the work array and trignometric tables used in transforming
!!  Also initialize the shapes that the arrays are going to take during
!!  the transform. In a three dimensional transform, the sequence is :
!!   Transform along X
!!   do a distributed transpose
!!   Transform along Y
!!    do a distributed transpose
!!  Transform along Z (in 2d, the last two steps are missing, but 
!!                     otherwise the process is similar)
!!  For all real transforms, the arrays go through three
!!  different shapes, the input shape is NX,NY/PY,NZ/PZ and is stored 
!!  in pfft_inLen. At the end of first distributed transpose, the shape is
!!  NY,NZ/PZ,NX/PY and is stored in pfft_midLen. After the second distributed
!!  transpose the shape is NZ,NX/PY,NY/PZ, and is stored in pfft_outLen.
!!  
!! ARGUMENTS
!!
!!   ndim - dimensionality of the problem
!!
!! NOTES
!!
!! DEV: My current understanding of how these things should be used: - KW   
!!
!!  WHEN / TYPE OF TRANFORM               DATA IS DESCRIBED BY
!! ======================================================================
!!
!!       |------------------------------- pfft_inLen,baseDatType(0)
!!       |
!! (first FFT,transformType(IAXIS))
!!       |
!!       |------------------------------- pfft_t1Len,baseDatType(IAXIS)
!!       |
!!  first Transpose
!!       |
!!       |------------------------------- pfft_midLen,baseDatType(IAXIS)
!!       |
!! (second FFT,transformType(JAXIS))
!!       |
!!       |------------------------------- pfft_t2Len,baseDatType(JAXIS)
!!       |
!!  second Transpose
!!       |
!!       |------------------------------- pfft_outLen,baseDatType(JAXIS)
!!       |
!! (third FFT,transformType(KAXIS))
!!       |
!!       |------------------------------- pfft_outLen,baseDatType(KAXIS)
!!       
!!
!! Currently we should always have baseDatType(KAXIS)==baseDatType(JAXIS).
!!
!! pfft_inlen              - IAXIS component measures lengths in terms of baseDatType.
!! pfft_midlen, pfft_outLen - IAXIS component measures lengths in terms of baseDatType.
!! pfft_t1len, pfft_t2Len - IAXIS component measures lengths in terms of baseDatType.
!! pfft_localLimits - 0-based, in terms of baseDatType (DEV: describes OUTPUT? of transform?!)
!! pfft_globalLen -            in terms of baseDatType (describes INPUT to the transform)
!!
!! wavesize (below)    <--  pfft_outLen, in terms of baseDatType.
!!***

subroutine gr_pfftInitMetaData(ndim)
#include "Pfft.h"
#include "constants.h"
  use gr_pfftData, ONLY : pfft_trigIaxis,pfft_trigJaxis, pfft_trigKaxis, &
       pfft_inLen,pfft_midLen, pfft_outLen, pfft_globalLen,&
       pfft_work1,pfft_work2, pfft_workSize, &
       pfft_ndim,pfft_comm,pfft_localLimits,&
       pfft_me,pfft_transformType,pfft_procGrid,pfft_scale,&
       pfft_t1Len,pfft_t2Len,pfft_wave,pfft_dimOrder, pfft_globalLen, pfft_myPE
  use gr_pfftInterface, ONLY : gr_pfftSetupDim,gr_pfftWave
  use Grid_interface, ONLY : Grid_getBlkIndexLimits
  use Driver_interface, ONLY : Driver_abortFlash

  implicit none

  integer, intent(IN) :: ndim
  integer,dimension(LOW:HIGH,MDIM) :: blkLimits,blkLimitsGC
  integer :: workArraySize, tempSize, i, blkFactor
  integer,dimension(MDIM) :: factor, blkLen
  integer :: len,waveSize,tempFact
  integer, dimension(MDIM) :: globalLen
  integer,dimension(MDIM) :: tfInDatSize, tfOutDatSizeDummy, datInc
  logical,save :: firstTime = .TRUE.

!!$  interface
!!$     integer function getMaxSize(inArray, outArray, factor, numProcs)
!!$       implicit none
!!$       integer, dimension(1:MDIM),intent(IN) :: inArray, outArray, factor
!!$       integer, intent(IN) :: numProcs
!!$     end function getMaxSize
!!$  end interface


  globalLen = pfft_globalLen

  call Grid_getBlkIndexLimits(1,blkLimits,blkLimitsGC)
  blkLen=blkLimits(HIGH,1:MDIM)-blkLimits(LOW,1:MDIM)+1


  !! Now allocate all the workspace needed
  !call gr_getComm(pfft_procGrid,pfft_comm,pfft_me) !NOT YET DEFINED
  i=IAXIS
  nullify(pfft_trigIaxis)
  factor(1:MDIM)=1


  !Allocating an array inside a subroutine and passing it back:
  !"http://www.sdsc.edu/~tkaiser/f90.html"
  call gr_pfftSetupDim(globalLen(i),pfft_transformType(i),&
       pfft_trigIaxis,pfft_scale(i),&
       tfInDatSize(i), tfOutDatSizeDummy(i),datInc(i),factor(i))

  pfft_dimOrder=1
  pfft_localLimits(:,:)=1
  pfft_midLen=1
  pfft_outLen=1

  if(ndim==1) then
     pfft_midLen = pfft_inLen
     pfft_outLen = pfft_inLen
     pfft_dimOrder(IAXIS)=IAXIS
  else
     i=JAXIS
     nullify(pfft_trigJaxis)
     call gr_pfftSetupDim(globalLen(i),pfft_transformType(i),&
          pfft_trigJaxis,pfft_scale(i),&
          tfInDatSize(i), tfOutDatSizeDummy(i),datInc(i),factor(i))
     if(ndim>2) then
        i=KAXIS
        nullify(pfft_trigKaxis)
        call gr_pfftSetupDim(globalLen(i),pfft_transformType(i),&
             pfft_trigKaxis,pfft_scale(i),&
             tfInDatSize(i), tfOutDatSizeDummy(i),datInc(i),factor(i))
     end if
  end if

  !! here find the largest size of work array
  !! First just find the size of the array at input
  if(ndim > 1) then

     !! This is the minimum required size for the work array
     workArraySize=pfft_inLen(IAXIS)*pfft_inLen(JAXIS)*pfft_inLen(KAXIS)

     !! Now find the maximum possible size after the first transpose.
     !! First step is to determine the sizes along the three dimensions after one
     !! transpose which changes the dimensions ordering from 1-2-3 to 2-3-1.

     !! along X axis, factor is 2 if the transform is complex.
     !! The assignment below is valid for 2D and 3D data.
     pfft_t1Len=pfft_inLen
     tempFact=max(tfInDatSize(JAXIS)/tfInDatSize(IAXIS),factor(IAXIS))

     !Chris: We have to consider whether the number of grid points in the IAXIS
     !is divisible by tempFact.
     !-------------------------------------------------------------
     if (mod(pfft_t1Len(IAXIS),tempFact) /= 0) then
        if (pfft_myPE == 0) then
           print *, "[gr_pfftInitMetadata]: Problem will fail... aborting."
           call Driver_abortFlash("[gr_pfftInitMetadata]: IAXIS not multiple of tempfact.")
        end if
     end if
     !-------------------------------------------------------------
     pfft_t1Len(IAXIS)=pfft_t1Len(IAXIS)/tempFact + datInc(IAXIS)
     pfft_midLen(IAXIS) = globalLen(JAXIS)

     if(ndim==2) then

        !! if ndim is 2, we only need to determine what is the local chunk when
        !! x - axis data is split among processors. We want this number to be
        !! Len/nproc if Len is a multiple of nproc, otherwise it should be
        !! Len/nproc + 1. Doing an integer operation (Len+nproc-1)/nproc results in
        !! the right answer.
        pfft_midLen(JAXIS) = (pfft_t1Len(IAXIS)+pfft_procGrid(JAXIS)-1)/pfft_procGrid(JAXIS)

        !! And in 2D, the final array sizes are the same as the sizes after one transpose
        pfft_outLen=pfft_midLen
        pfft_dimOrder(IAXIS)=JAXIS
        pfft_dimOrder(JAXIS)=IAXIS
     else !! but if the data is 3 dimensional

        !! Since the transformation is from Nx,Ny/Py,Nz/Pz to Ny,Nz/Pz,Nx/Py
        !! the parallel distribution of the K axis doesn't change, it only 
        !! becomes the second axis of the array
        pfft_midLen(JAXIS)=pfft_t1Len(KAXIS)

        !! What was IAXIS, now becomes the third axis which is divided among Py processors
        pfft_midLen(KAXIS) = (pfft_t1len(IAXIS)+&
             pfft_procGrid(JAXIS)-1)/pfft_procGrid(JAXIS)


        !! The second transpose step applies only in 3D cases.
        !! In this transpose step, the data distribution changes from Ny,Nz/Pz,Nx/Py to
        !! Nz,Nx/Py,Ny/Pz. Intermediate values calculated in the midLen step are used
        !! The logic followed is similar to that in the first transpose using midLen
        !! values as input
        pfft_t2Len=pfft_midLen
        tempFact=max(tfInDatSize(KAXIS)/tfInDatSize(JAXIS),factor(JAXIS))
        pfft_t2Len(IAXIS)=pfft_t2Len(IAXIS)/tempFact + datInc(JAXIS)

        pfft_outLen(IAXIS)=globalLen(KAXIS)
        pfft_outLen(JAXIS)=pfft_midLen(KAXIS)
        pfft_outLen(KAXIS)=(pfft_t2Len(IAXIS)+pfft_procGrid(KAXIS)-1)/pfft_procGrid(KAXIS)
        pfft_dimOrder(IAXIS)=KAXIS
        pfft_dimOrder(JAXIS)=IAXIS
        pfft_dimOrder(KAXIS)=JAXIS

     end if

     !! Now that all possible configurations of the transforms sizes are know, we can
     !! calculate the work array size
     factor(IAXIS) = max(factor(IAXIS),tfInDatSize(IAXIS))

     tempSize = pfft_midLen(IAXIS)*pfft_midLen(KAXIS)*pfft_midLen(JAXIS)
     tempSize=tempSize*min(2,factor(IAXIS)*factor(JAXIS)*factor(KAXIS))
     workArraySize=max(workArraySize,tempSize)

     tempSize=pfft_outLen(IAXIS)*pfft_outLen(JAXIS)*pfft_outLen(KAXIS)
     tempSize=tempSize*min(2,factor(IAXIS)*factor(JAXIS)*factor(KAXIS))
     workArraySize=max(workArraySize,tempSize)

     !Chris: We need to consider the workArraySize from the transpose 
     !perspective.  This sizes the work arrays such that the transpose 
     !can proceed without the MPI_ALLTOALL being accessed out of bounds, 
     !and the unit test works (e.g. 2D and nx=24,ny=24,iProcs=3,jProcs=3).
     !-------------------------------------------------------------
     !This is the calculation that is performed in the transpose.
     tempSize = max(tempSize, getMaxSize(pfft_t1Len, pfft_midLen, factor, pfft_procGrid(JAXIS)))
     if (pfft_ndim == 3) then
        tempSize = max(tempSize, getMaxSize(pfft_t2Len, pfft_outLen, factor, pfft_procGrid(KAXIS)))
     end if
     if (tempsize > workArraySize) then
        if (pfft_myPE == 0) then
           print *, "[gr_pfftInitMetadata]: WARNING... making work arrays larger artificially!!!", &
                " Size was:", workArraySize, "now:", tempSize
        end if
        workArraySize=max(workArraySize,tempSize)
     end if
     !-------------------------------------------------------------

  else
     workArraySize=1
  end if
  pfft_workSize = workArraySize

  allocate(pfft_work1(workArraySize))
  allocate(pfft_work2(workArraySize))
  pfft_work1=0.0
  pfft_work2=0.0

  if (pfft_myPE == 0) then
     if (firstTime) then
999     format(5(1x,A12,3I12:/),1(1x,A14,I11:/))
        print 999, "pfft_inLen:", pfft_inLen, &
             "pfft_midLen:", pfft_midLen, &
             "pfft_outLen:", pfft_outLen, &
             "pfft_t1Len:", pfft_t1Len, &
             "pfft_t2Len:", pfft_t2Len, &
             "workarraysize:", workarraysize
     end if
     firstTime = .FALSE.
  end if


  waveSize=0
  do i = 1,NDIM
     len=pfft_outLen(i)
!!$     if((pfft_transformType(i)==PFFT_REAL2C).or.(pfft_transformType(i)==PFFT_REAL2C_STUFF).or.&
!!$          (pfft_transformType(i)==PFFT_REAL2C_EXTEND).or.&
!!$          (pfft_transformType(i)==PFFT_COMPLEX))len=len*2
     !waveSize=max(waveSize,pfft_outLen(i))
     waveSize=max(waveSize,len)
     call gr_pfftGetLocalLimits(i,pfft_dimOrder(i))
  end do
  allocate(pfft_wave(waveSize,MDIM))
  pfft_wave=0.0  !We must initialise here, as not all elements are written to in gr_pfftWave.
  call gr_pfftWave()

  return

contains

  integer function getMaxSize(inArray, outArray, factor, numProcs)

#include "Flash.h"

    implicit none
    integer, dimension(1:MDIM),intent(IN) :: inArray, outArray, factor
    integer, intent(IN) :: numProcs
    integer :: tempSize

    tempSize = inArray(JAXIS) * inArray(KAXIS) * outArray(NDIM) * numProcs
    tempSize = tempSize * min(2,factor(IAXIS)*factor(JAXIS)*factor(KAXIS))
    getMaxSize = tempSize
end function getMaxSize

end subroutine gr_pfftInitMetaData
