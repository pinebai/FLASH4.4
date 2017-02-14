!!****if* source/Grid/GridBoundaryConditions/OneRow/Flash3HSE/gr_applyFlash3HSEBC
!!
!! NAME
!!
!!  gr_applyFlash3HSEBC
!!
!! SYNOPSIS
!!
!!  call gr_applyFlash3HSEBC(integer(IN) :: bcType,
!!                           integer(IN) :: bcDir,
!!                           integer(IN) :: guard,
!!                           real(INOUT) :: dataRow(2*guard,NUNK_VARS),
!!                           integer(IN) :: face,
!!                           real(IN)    :: cellCenterSweepCoord(*),
!!                           real(IN)    :: secondCoord,
!!                           real(IN)    :: thirdCoord)
!!
!! DESCRIPTION
!!
!!   Implementation of hydrostatic boundary conditions according to Zingal, Dursi, et al,
!!   based on implementation by Dean Townsley.
!!
!! ARGUMENTS
!!
!!   bcType : the type of boundary condition being applied.
!!            Should be a negative number selecting one of the variants
!!            of HSE BCs.
!!
!!   bcDir :  the dimension along which to apply boundary conditions,
!!            can take values of IAXIS, JAXIS, and KAXIS.
!!
!!   guard :  number of guardcells
!!
!!   dataRow: storage for the data being operated upon.
!!            one half is input, the other half is output.
!!            Which is which depends on face.
!!
!!   face :   can take values LOW and HIGH, defined in constants.h,
!!            to indicate whether to apply boundary on lowerface or 
!!            upperface
!!
!!  cellCenterSweepCoord : vector of (at least) 2*guard cell center coordinate
!!                         values in the bcDir direction. The elements of this
!!                         array give the locations of corresponding data in
!!                         dataRow.
!!
!!  secondCoord,thirdCoord: scalar cell center coordinate values in the coordinate
!!                         directions perpendicular to the direction given by bcDir.
!!                         This is used in this implementation only for passing
!!                         the information along to Gravity_accelAtCoords, which
!!                         is the Gravity interface that is invoked to get the
!!                         strength of the gravitational field.
!!                         The meaning depends on the sweep direction bcDir as
!!                         follows
!!                          bcDir   |    secondCoord       thirdCoord
!!                          ------------------------------------------
!!                          IAXIS   |    Y(j) *            Z(k) **
!!                          JAXIS   |    X(i)              Z(k) **
!!                          KAXIS   |    X(i)              Y(j)
!!                         *)  if NDIM > 1
!!                         **) if NDIM > 2
!!                         These dummy arguments are ignored (and an implementation
!!                         of this interface should not attempt to access them) if
!!                         they do not make sense based on the dimensionality of
!!                         the problem.
!!
!! NOTES
!!
!!
!!  This implementation may or may not be appropriate for a given application.
!!
!! HISTORY
!!
!!  Dean Townsley
!!  Dongwook Lee
!!  Klaus Weide
!!
!!***


subroutine gr_applyFlash3HSEBC(bcType,bcDir,guard,dataRow,face,&
     cellCenterSweepCoord, secondCoord,thirdCoord)

!!$subroutine Grid_bcApplyToRegionSpecialized(bcType,gridDataStruct,&
!!$     guard,axis,face,region,regionSize,mask,applied,&
!!$     blockHandle,secondDir,thirdDir,endPoints,blkLimitsGC, idest)

#include "constants.h"
#include "Flash.h"
#include "Eos.h"

  use Grid_data, ONLY :gr_meshMe,gr_domainBC
  use Driver_interface, ONLY : Driver_abortFlash
  use gr_bcInterface, ONLY : gr_hseStep
  use gr_bcHseData, ONLY : gr_bcHseGravConst, HSE_FORWARD, HSE_BACKWARD, HSE_SETTEMP
  use Eos_interface, ONLY : Eos

  implicit none

  integer,intent(IN):: bcType,bcDir,guard,face
  real,dimension(2*guard,NUNK_VARS),intent(INOUT)::dataRow
  real,intent(IN):: cellCenterSweepCoord(*), secondCoord,thirdCoord
!!$  integer, intent(IN) :: bcType,axis,face,guard,gridDataStruct
!!$  integer,dimension(REGION_DIM),intent(IN) :: regionSize
!!$  real,dimension(regionSize(BC_DIR),&
!!$       regionSize(SECOND_DIR),&
!!$       regionSize(THIRD_DIR),&
!!$       regionSize(STRUCTSIZE)),intent(INOUT)::region
!!$  logical,intent(IN),dimension(regionSize(STRUCTSIZE)):: mask
!!$  logical, intent(OUT) :: applied
!!$  integer,intent(IN) :: blockHandle
!!$  integer,intent(IN) :: secondDir,thirdDir
!!$  integer,intent(IN),dimension(LOW:HIGH,MDIM) :: endPoints, blkLimitsGC
!!$  integer,intent(IN),OPTIONAL:: idest

  integer :: i,j,k, sizeGC, istart, iend, step, direction
  integer :: velVarBcDir, velVarSecondDir, velVarThirdDir
  real    :: deltax
  real, dimension(EOS_NUM) :: eosData


!=====================================================================

#ifndef FLASH_EOS
  call Driver_abortFlash('Cannot execute gr_applyFlash2HSEBC without the Eos unit!')
#endif

  select case (bcDir)
  case(IAXIS)
     velVarBcDir =     VELX_VAR
     velVarSecondDir = VELY_VAR
     velVarThirdDir =  VELZ_VAR
  case(JAXIS)
     velVarBcDir =     VELY_VAR
     velVarSecondDir = VELX_VAR
     velVarThirdDir =  VELZ_VAR
  case(KAXIS)
     velVarBcDir =     VELZ_VAR
     velVarSecondDir = VELX_VAR
     velVarThirdDir =  VELY_VAR
  end select

  ! We are going to assume that gravity is along the bcDir direction.
  ! apply our hydrostatic in both directions


  ! assume this is uniform
  deltax = cellCenterSweepCoord(2)-cellCenterSweepCoord(1)


        !-------------------
        !  first take care of the velocities
        if (face==HIGH) then
           do i = 1,guard
              ! zero-gradient everything (cannot know what the user has defined that needs to be propagated)
              dataRow(guard+i,:) = dataRow(guard,:)
              if (bcType == HYDROSTATIC_NVREFL) then
                 dataRow(guard+i,velVarBcDir) = -dataRow(guard+1-i,velVarBcDir)
                 dataRow(guard+i,velVarSecondDir) =  dataRow(guard+1-i,velVarSecondDir)
                 dataRow(guard+i,velVarThirdDir) = dataRow(guard+1-i,velVarThirdDir)
              else if (bcType == HYDROSTATIC_NVOUT) then
                 dataRow(guard+i,velVarBcDir) = dataRow(guard,velVarBcDir)
                 dataRow(guard+i,velVarSecondDir) = dataRow(guard,velVarSecondDir)
                 dataRow(guard+i,velVarThirdDir) = dataRow(guard,velVarThirdDir)
              else if (bcType == HYDROSTATIC_NVDIODE) then
                 dataRow(guard+i,velVarBcDir) = max(0.0,dataRow(guard,velVarBcDir))
                 dataRow(guard+i,velVarSecondDir) = dataRow(guard,velVarSecondDir)
                 dataRow(guard+i,velVarThirdDir) = dataRow(guard,velVarThirdDir)
              endif
          enddo
        else if (face==LOW) then
           do i = 1,guard
              ! zero-gradient everything (cannot know what the user has defined that needs to be propagated)
              dataRow(i,:)     = dataRow(guard+1,:)
              if (bcType == HYDROSTATIC_NVREFL) then
                 dataRow(i,velVarBcDir) = -dataRow(2*guard+1-i,velVarBcDir)
                 dataRow(i,velVarSecondDir) =  dataRow(2*guard+1-i,velVarSecondDir)
                 dataRow(i,velVarThirdDir) =  dataRow(2*guard+1-i,velVarThirdDir)
              else if (bcType == HYDROSTATIC_NVOUT) then
                 dataRow(i,velVarBcDir) = dataRow(guard+1,velVarBcDir)
                 dataRow(i,velVarSecondDir) = dataRow(guard+1,velVarSecondDir)
                 dataRow(i,velVarThirdDir) = dataRow(guard+1,velVarThirdDir)
              else if (bcType == HYDROSTATIC_NVDIODE) then
                 dataRow(i,velVarBcDir) = min(0.0,dataRow(guard+1,velVarBcDir))
                 dataRow(i,velVarSecondDir) = dataRow(guard+1,velVarSecondDir)
                 dataRow(i,velVarThirdDir) = dataRow(guard+1,velVarThirdDir)
              endif
          enddo
        endif
       
        if (face==HIGH) then
           ! fill stuff out
           istart = guard+1
           iend   = 2*guard
           step   = 1
           direction = HSE_FORWARD
        else
           istart = guard
           iend   = 1
           step   = -1
           direction = HSE_BACKWARD
        endif
        do i = istart, iend, step
           ! density guess value
           dataRow(i,DENS_VAR) = dataRow(i-step,DENS_VAR)
           call gr_hseStep(dataRow(:,DENS_VAR), &
                             dataRow(:,TEMP_VAR), &
                             dataRow(:,YE_MSCALAR),   &
                             dataRow(:,SUMY_MSCALAR), &
                             i,gr_bcHseGravConst, deltax, direction,2, HSE_SETTEMP)
      
           ! now get all the eos stuff and fill it in
           eosData(EOS_DENS) = dataRow(i,DENS_VAR)
           eosData(EOS_TEMP) = dataRow(i,TEMP_VAR)
           eosData(EOS_ABAR) = 1.0/dataRow(i,SUMY_MSCALAR)
           eosData(EOS_ZBAR) = dataRow(i,YE_MSCALAR)*eosData(EOS_ABAR)
           call Eos(MODE_DENS_TEMP, 1, eosData)
           dataRow(i,PRES_VAR) = eosData(EOS_PRES)
           dataRow(i,EINT_VAR) = eosData(EOS_EINT)
           dataRow(i,GAME_VAR) = eosData(EOS_PRES)/(eosData(EOS_EINT)*eosData(EOS_DENS)) +1.0
           dataRow(i,GAMC_VAR) = eosData(EOS_GAMC)

           ! Dongwook: put 1/2 factor for energy
           dataRow(i,ENER_VAR) = eosData(EOS_EINT) + 0.5*(dataRow(i,VELX_VAR)**2 &
                                      + dataRow(i,VELY_VAR)**2 + dataRow(i,VELZ_VAR)**2)

        enddo

  return
end subroutine gr_applyFlash3HSEBC
