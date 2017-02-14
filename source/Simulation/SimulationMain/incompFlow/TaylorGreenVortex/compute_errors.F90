!!****if* source/Simulation/SimulationMain/incompFlow/TaylorGreenVortex/compute_errors
!!
!! NAME
!!
!!  compute_errors
!!
!! SYNOPSIS
!!
!!  call compute_errors(:: mype,
!!                      integer, INTENT(IN)  :: blockcount,
!!                      integer, INTENT(IN), dimension(blockCount)  :: blocklist)
!!
!! DESCRIPTION
!!
!!
!! ARGUMENTS
!!
!!   mype : my Processor Number
!!
!!   blockcount : 
!!
!!   blocklist : 
!!
!! AUTOGENROBODOC
!!
!!
!!***

! subroutine compute_errors
! Computes velocity and pressure error for Taylor Green Vortex problem.
! --------------------------------------------------------------------------

  subroutine compute_errors(mype,blockCount,blockList)

#include "Flash.h"

  ! Modules Use:
  use Grid_interface, ONLY : Grid_getDeltas, Grid_getBlkBC,       &
      Grid_getBlkPtr, Grid_releaseBlkPtr, Grid_getBlkIndexLimits, &
      Grid_getBlkBoundBox, Grid_getBlkCenterCoords
  use gr_interface, ONLY : gr_findMean
  use Driver_data, ONLY : dr_simTime,dr_nstep
  use IncompNS_data, ONLY : ins_invRe
  use Simulation_data, ONLY : sim_uconv, sim_vconv
      
#ifdef  FLASH_GRID_PARAMESH
  use tree, only : lrefine,lrefine_max
#endif
  

  implicit none
#include "constants.h"
#include "IncompNS.h"
  include 'mpif.h'
  ! Local Variables:
  integer :: mype
  integer, INTENT(IN) :: blockCount
  integer, INTENT(IN), dimension(blockCount) :: blockList

  integer, parameter ::  ng = NGUARD
  integer, parameter ::  nxi= NGUARD + NXB
  integer, parameter ::  nyj= NGUARD + NYB
  integer, parameter ::  nzk= NGUARD + NZB
  integer, parameter ::  nxc= NGUARD + NXB + 1
  integer, parameter ::  nyc= NGUARD + NYB + 1
  integer, parameter ::  nzc= NGUARD + NZB + 1

  real, pointer, dimension(:,:,:,:) :: solnData, facexData,faceyData
  
  real del(3),dx,dy,dz

  real meanPres

  integer lb,blockID,ierr

  real erru_inf,errv_inf,errp_inf,erru_L2,errv_L2,errp_L2
  integer Nux,Nvy,Npp
  real    :: mvisc

  real xedge,xcell,yedge,ycell

  integer i,j

  real coord(MDIM),bsize(MDIM)
  real, dimension(2,MDIM) :: boundBox

  logical, save :: firstime = .TRUE.

  ! -------------------------------------------------------------------

  mvisc = ins_invRe

  erru_inf =0.
  erru_L2  =0.
  Nux = 0
  errv_inf =0.
  errv_L2  =0.
  Nvy = 0
  errp_inf =0.
  errp_L2  =0.
  Npp = 0

  call gr_findMean(PRES_VAR,2,.false.,meanPres)

  write(*,*) 'Mean Pressure =',meanPres
  do lb = 1,blockCount

     blockID = blockList(lb)

#ifdef  FLASH_GRID_PARAMESH
     if (lrefine(blockID) .eq. lrefine_max) then
#endif

     ! Get blocks dx, dy ,dz:
     call Grid_getDeltas(blockID,del)

     ! Get Coord and Bsize for the block:
     ! Bounding box:
     call Grid_getBlkBoundBox(blockId,boundBox)
     bsize(:) = boundBox(2,:) - boundBox(1,:)

     call Grid_getBlkCenterCoords(blockId,coord)

     ! Point to blocks center and face vars:
     call Grid_getBlkPtr(blockID,facexData,FACEX)
     call Grid_getBlkPtr(blockID,faceyData,FACEY)
     call Grid_getBlkPtr(blockID,solnData,CENTER)


     do j = NGUARD,NYB+2*NGUARD
        yedge = coord(JAXIS) - bsize(JAXIS)/2.0 +    &
                real(j - NGUARD - 1)*del(JAXIS) - sim_vconv*dr_simTime
        ycell = yedge + del(JAXIS)/2.0

        do i = NGUARD,NXB+2*NGUARD
           xedge = coord(IAXIS) - bsize(IAXIS)/2.0 + &
                   real(i - NGUARD - 1)*del(IAXIS) - sim_uconv*dr_simTime
           xcell = xedge + del(IAXIS)/2.0

           facexData(EVEL_FACE_VAR,i,j,1) = -EXP(-2.0*mvisc*dr_simTime)* &
                                             COS(xedge)*SIN(ycell)

           faceyData(EVEL_FACE_VAR,i,j,1) =  EXP(-2.0*mvisc*dr_simTime)* &
                                             SIN(xcell)*COS(yedge)

           solnData(EPRS_VAR,i,j,1)  = -0.25*EXP(-4.0*mvisc*dr_simTime)* &
                        ( COS(2.0*xcell) + COS(2.0*ycell) ) + meanPres 

        enddo
     enddo
     


     erru_L2  = erru_L2 + &
                     (sum((facexData(VELC_FACE_VAR,ng+1:nxc,ng+1:nyj,1) - sim_uconv - &
                           facexData(EVEL_FACE_VAR,ng+1:nxc,ng+1:nyj,1))**2.))
     Nux = Nux + (NXB+1)*NYB


     erru_inf = max(erru_inf, &
                maxval(abs(facexData(VELC_FACE_VAR,ng+1:nxc,ng+1:nyj,1) - sim_uconv - &
                           facexData(EVEL_FACE_VAR,ng+1:nxc,ng+1:nyj,1))))


     errv_L2  = errv_L2 + &
                     (sum((faceyData(VELC_FACE_VAR,ng+1:nxi,ng+1:nyc,1) - sim_vconv - &
                           faceyData(EVEL_FACE_VAR,ng+1:nxi,ng+1:nyc,1))**2.))
     Nvy = Nvy + NXB*(NYB+1)


     errv_inf = max(errv_inf, &
                maxval(abs(faceyData(VELC_FACE_VAR,ng+1:nxi,ng+1:nyc,1) - sim_vconv - &
                           faceyData(EVEL_FACE_VAR,ng+1:nxi,ng+1:nyc,1))))


     errp_L2  = errp_L2 + &
                     (sum((solnData(PRES_VAR,ng+1:nxi,ng+1:nyj,1)- &
                           solnData(EPRS_VAR,ng+1:nxi,ng+1:nyj,1))**2.))
     Npp = Npp + NXB*NYB


     errp_inf = max(errp_inf, &
                maxval(abs(solnData(PRES_VAR,ng+1:nxi,ng+1:nyj,1)- &
                           solnData(EPRS_VAR,ng+1:nxi,ng+1:nyj,1))))



     ! Release pointers:
     call Grid_releaseBlkPtr(blockID,facexData,FACEX)
     call Grid_releaseBlkPtr(blockID,faceyData,FACEY)
     call Grid_releaseBlkPtr(blockID,solnData,CENTER)
   

#ifdef  FLASH_GRID_PARAMESH
     endif
#endif

        
  enddo

  ! Here MPI_allreduce sum computed errors:
  



  write(*,*) 'Nux,Nvy,Npp=',Nux,Nvy,Npp

  if (Nux.NE.0) erru_L2  = sqrt(erru_L2/real(Nux))
  if (Nvy.NE.0) errv_L2  = sqrt(errv_L2/real(Nvy))
  if (Npp.NE.0) errp_L2  = sqrt(errp_L2/real(Npp))

  ! Gather error information:



  if (mype .eq. 0) then
  !write(*,*) 'U Velocities:'
  write(*,*) mype,'einf U=',erru_inf,'eL2 U=',erru_L2
  !write(*,*) 'V Velocities:'
  write(*,*) mype,'einf V=',errv_inf,'eL2 U=',errv_L2
  !write(*,*) 'P Presures:'
  write(*,*) mype,'einf P=',errp_inf,'eL2 P=',errp_L2
  endif

  ! Write error files:
  if (firstime) then

 !    open(unit=33, file='./IOData/uerror_time.res', status='replace')
 !    open(unit=34, file='./IOData/verror_time.res', status='replace')
 !    open(unit=35, file='./IOData/perror_time.res', status='replace')    

     firstime = .false. 

  else

 !    open(unit=33, file='./IOData/uerror_time.res', status='old', position='append')
 !    open(unit=34, file='./IOData/verror_time.res', status='old', position='append')
 !    open(unit=35, file='./IOData/perror_time.res', status='old', position='append')   

  endif

 ! write(33,*) dr_nstep,erru_inf,erru_L2
 ! write(34,*) dr_nstep,errv_inf,errv_L2
 ! write(35,*) dr_nstep,errp_inf,errp_L2

 ! close(33)
 ! close(34)
 ! close(35) 


  return
  end subroutine


