!!****if* source/Grid/GridSolvers/Pfft/gr_pfftSetupDim
!!
!! NAME
!!
!!  gr_pfftSetupDim
!!
!! SYNOPSIS
!!  
!!  call gr_pfftSetupDim(integer(IN) :: length,
!!                       integer(IN) :: transformType,
!!                       real(OUT)   :: trig(:),
!!                       real(OUT)   :: scale,
!!                       integer(OUT):: inDatSize,
!!                       integer(OUT):: outDatSize,
!!                       integer(OUT):: datInc,
!!                       integer(OUT):: factor)
!! 
!! DESCRIPTION
!!
!!   This routine computes the trignometric tables and scaling
!!   factors for each individual dimension in the transform. 
!!
!! ARGUMENTS
!!  
!!  length - length of vectors which will be transformed
!!  transformType - the allowed value for this argument are:
!!                 PFFT_REAL2C  (real to complex)
!!                 PFFT_REAL  (real to real)
!!                 PFFT_COMPLEX (complex to compex)
!!                 PFFT_COS (cosine)
!!                 PFFT_SIN (sine)
!!                 PFFT_COSQ (cos with a phase shift)
!!                 PFFT_SINQ (sin with a phase shift)
!!  trig  - this array contains precomputed trignometric tables
!!  scale - the scaling factor
!!  inDatSize - this argument returns a value 2 if the input data of
!!           transform is interpreted as complex,
!!           otherwise the value is 1.
!!  outDatSize - this argument returns a value 2 if the output data of
!!           transform is interpretable as complex,
!!           otherwise the value is 1.
!!  datInc - number of additional elements for which space is required
!!           in the output array, in multiples of outDatSize; 0 in
!!           most cases.
!!  factor - this argument returns a value 2 if transform output is
!!           complex (i.e., outDatSize=2) while input is real data;
!!           otherwise the value is 1.
!!
!!***

subroutine gr_pfftSetupDim(length, transformType, trig, scale, inDatSize, outDatSize, datInc, factor)
#include "Pfft.h"
  implicit none
  integer,intent(IN) :: length,transformType
  real,dimension(:), pointer :: trig 
  real,intent(OUT) :: scale
  integer,intent(OUT) :: inDatSize, outDatSize
  integer,intent(OUT) :: datInc, factor
  integer :: trigLen

  if(associated(trig)) deallocate(trig)
  nullify(trig)

  inDatSize=1
  outDatSize=1
  factor=1
  datInc=0
  select case(transformType)
  case(PFFT_COSQ,PFFT_COS_CC)
     scale = 0.25/length
     trigLen=3*length+15
     allocate(trig(trigLen))
     call cosqi(length,trig)
  case(PFFT_SINQ,PFFT_SIN_CC)
     scale = 0.25/length
     trigLen=3*length+15
     allocate(trig(trigLen))
     call sinqi(length,trig)
  case(PFFT_COS_IV)
     scale = 0.03125/length        !DEV: ?
     trigLen=6*length+15
     allocate(trig(trigLen))
     call cosqi(2*length,trig)
  case(PFFT_SIN_IV)
     scale = 0.03125/length        !DEV: ?
     trigLen=6*length+15
     allocate(trig(trigLen))
     call sinqi(2*length,trig)
  case(PFFT_SIN)
     scale = 0.5/(length+1)
     trigLen=(5*length+30)/2
     allocate(trig(trigLen))
     !call sinti(length,trig)  !f77 fftpack
     call rsinti(length,trig)  !f90 fftpack
  case(PFFT_COS)
     scale = 0.5/(length-1)
     trigLen=3*length+15
     allocate(trig(trigLen))
     !call costi(length,trig)  !f77 fftpack
     call rcosti(length,trig)  !f90 fftpack
  case(PFFT_REAL,PFFT_REAL2C,PFFT_REAL2C_STUFF,PFFT_REAL2C_EXTEND)
     scale = 1.0/length
     trigLen=2*length+15
     allocate(trig(trigLen))
     call rffti(length,trig)
     if(transformType .NE. PFFT_REAL) then
        outDatSize=2
        factor=2
        if(transformtype==PFFT_REAL2C_EXTEND) datInc=1
     end if
  case(PFFT_COMPLEX)
     scale = 1.0/length
     trigLen=4*length+15
     allocate(trig(trigLen))
     call cffti(length,trig)
     inDatSize=2
     outDatSize=2
  end select
  return
end subroutine gr_pfftSetupDim
