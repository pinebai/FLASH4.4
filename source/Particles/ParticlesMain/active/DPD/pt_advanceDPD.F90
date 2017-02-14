!!****if* source/Particles/ParticlesMain/active/DPD/pt_advanceDPD
!!
!! NAME
!!
!!  pt_advanceDPD
!!
!! SYNOPSIS
!!
!!  pt_advanceDPD(real(IN) :: dtOld,
!!                  real(IN) :: dtNew,
!!                  real(inout):: particles(:,p_count),
!!                  integer(in):: p_count)
!!                  integer(in):: ind)
!!
!! ARGUMENTS
!!  
!!   dtOld : previous time interval 
!!   dtNew : current time interval
!!   particles -- particles to advance
!!   p_count  -- the number of particles in the list to advance
!!   ind    --- index for type into pt_typeInfo
!!  
!! DESCRIPTION
!!
!!  Time advancement routine for the particle module.
!!  
!!  This version is the forward Euler advancement for the active 
!!  submodule.
!!
!!***

!===============================================================================

subroutine pt_advanceDPD (dtOld,dtNew,particles,p_count, ind)
    
  use Particles_data, ONLY: useParticles,pt_meshMe,pt_indexList, pt_indexCount, pt_maxperproc, pt_numlocal
  use Timers_interface, ONLY : Timers_start, Timers_stop
  use Driver_interface, ONLY : Driver_getNStep
  use pt_dpdData, ONLY : pt_dpdUpdateCycle, pt_dpdNstep,pt_dpdLambda

  implicit none

#include "Flash.h"
#include "constants.h"
    
  real, INTENT(in) :: dtOld, dtNew
  integer,intent(in) :: p_count, ind
  real,dimension(NPART_PROPS,p_count),intent(inout) :: particles

!-------------------------------------------
  logical :: regrid = .false.
  integer :: i,j
  integer,dimension(MDIM):: opind,npind,ovind,intvind,nvind,ofind,nfind
  integer :: particlesPerBlk
  real,dimension(NDIM,p_count)::pos,npos,v,vnew,f,fnew,v_telda
  real :: dt
  real,dimension(NPART_PROPS,p_count) :: particles2
 
!-------------------------------------------  

  if (.not.useParticles ) return

  call Timers_start ("particles")
  call Driver_getNStep(pt_dpdNstep)
  
  dt=dtOld;
  particles2=particles;
    
  !write(*,*)pt_dpdupdateCycle,pt_dpdNstep
  ! This is useful if the particles are not to be moved every time step.
  if((pt_dpdupdateCycle==1).or.(MOD(pt_dpdNstep,pt_dpdupdateCycle)==0)) then     
     call Grid_moveParticles(particles,NPART_PROPS,pt_maxPerProc,pt_numLocal, &
          pt_indexList, pt_indexCount, regrid) 
  end if
  
  ! Now update the pt_typeInfo data structure
  call pt_updateTypeDS(particlesPerBlk)
  
  ! Set the indices for variables at n and n+1 
  call pt_dpdSetIndices(opind,npind,ovind,nvind,intvind,ofind,nfind)
  
  ! Initialize the fnew,vnew,
  fnew(:,:)=0;
  vnew(:,:)=0;
  
  ! Read the data from the particles array to update the positions
  
  f(1:NDIM,1:p_count)      = particles(nfind(1:NDIM),1:p_count);
  v(1:NDIM,1:p_count)      = particles(nvind(1:NDIM),1:p_count);
  pos(1:NDIM,1:p_count)    = particles(npind(1:NDIM),1:p_count);
  v_telda(1:NDIM,1:p_count)= particles(intvind(1:NDIM),1:p_count);
  
  
  ! Update forces: rest of the previous step  
  call pt_dpdNonBondedForces(pos,v_telda,particles(BDT_PART_PROP,1:p_count), & 
       particles(PRNT_PART_PROP,1:p_count),&
       particles(BDYT_PART_PROP,1:p_count), &
       particles(INTR_PART_PROP,1:p_count),fnew,p_count);
  
  !Update velocity 1st time
  vnew= v + 0.5 *  dt * (f+fnew);
  
  ! The first half of inegration for the new step 
  npos=pos+ dt*vnew + 0.5 * dt*dt * fnew ;
  v_telda= vnew + pt_dpdLambda * dt * fnew;
  
     ! Update the particles array
!!$  write(*,*)opind(1:NDIM),POSX_PART_PROP,POSZ_PART_PROP
!!$  stop
     particles(npind(1:NDIM),1:p_count)= npos(1:NDIM,1:p_count);
     particles(opind(1:NDIM),:)= pos(1:NDIM,1:p_count);
     particles(ofind(1:NDIM),:)= f(1:NDIM,1:p_count);
     particles(nfind(1:NDIM),:)= fnew(1:NDIM,1:p_count);
     particles(intvind(1:NDIM),:)= v_telda(1:NDIM,1:p_count);
     particles(nvind(1:NDIM),:)= vnew(1:NDIM,1:p_count);
     particles(ovind(1:NDIM),:)= v(1:NDIM,1:p_count);

  
  
  
#ifdef DEBUG_VPARTICLES
  if (pt_meshMe==3) then
     write(*,*)
     do i=1,p_count
        write(*,*)'Particle',i,'X=',particles(POSX_PART_PROP,i), &
             'Y=',particles(POSY_PART_PROP,i), &
             'Z=',particles(POSZ_PART_PROP,i), &
             'TAG=',particles(TAG_PART_PROP,i)
     end do
  end if
#endif
  
  ! Remove the virtual particles after the forcing
  ! Local number of points
  pt_numlocal=p_count
  call Particles_clean()
  
#ifdef DEBUG_VPARTICLES2
  !if (pt_globalMe==3) then
  write(*,*)'pt_numlocal',pt_numlocal,'pt_globalMe',pt_globalMe
  do i=1,pt_numlocal
     write(*,*)'Particle',i,'X=',particles(POSX_PART_PROP,i), &
          'Y=',particles(POSY_PART_PROP,i), &
          'Z=',particles(POSZ_PART_PROP,i), &
          'TAG=',particles(TAG_PART_PROP,i),&
          'Proc=',pt_globalMe
  end do
  !end if
#endif
  
  call Timers_stop ("particles")
  
  return
  
end subroutine pt_advanceDPD

!===============================================================================

