!!****if* source/Particles/ParticlesMain/active/DPD/pt_dpdNonBondedForces
!!
!! NAME
!!    pt_DPDupdateForces
!!
!! SYNOPSIS
!!    pt_dpdNonBondedForces( real,   INTENT(IN)     ::pos
!!                         real,   INTENT(IN)     ::v
!!                         real,   INTENT(IN)     ::btypes
!!                         real,   INTENT(IN)     ::parent 
!!                         real,   INTENT(IN)     ::parentType
!!                         real  , INTENT(OUT)    ::fvec
!!                         integer,INTENT(IN)     ::p_count)
!!
!! DESCRIPTION
!!
!!    Calculates the bonded and nonbonded forces between particles on each processor
!!    The Aij array holds the repulsion parameters between all beads types
!!    The property BDT_PART_PROP for each particle contains its bead type, 
!!    therefore, for calculating the force between the two particles i and j  
!!    The value stored in Aij(particle(BDT_PART_PROP,i),particle(BDT_PART_PROP,j)) will be used  
!!    to evaluate the forces between these two particles
!!
!!
!!
!! ARGUMENTS
!!
!!  particles :A real type array holding the local particles (real and virtual) on this processor
!! 
!!            for this version of the routine
!!
!! updateRefine : is true if the routine wished to retain the already
!!                initialized particles instead of reinitializing them
!!                as the grid refine.
!!
!!***

subroutine pt_dpdNonBondedForces(pos,v,btypes,parents,parentType, &
     internalIndex,fvec,p_count)
  
  use Simulation_data
  use Driver_data,ONLY:dr_globalMe,dr_nstep
  implicit none
#include "Flash.h"
#include "constants.h"

  integer,INTENT(IN) :: p_count
  real,dimension(NDIM,p_count),INTENT(IN) :: pos,v
  real,dimension(p_count),INTENT(IN) :: btypes,parents,parentType,internalIndex
  real,dimension(NDIM,p_count),INTENT(OUT) :: fvec

  !---- Local variables
  real :: vij,l,fval,a_rep,r_dot_v
  real :: lij(NDIM),vij3(NDIM)
  real :: zeta,Fij_C,Fij_R,Fij_D
  integer :: np,i,j,BodyType
  integer :: Local_numBodies,nLinks,bead1,bead2
  logical :: Res(MAXCONN)
  !-------------------------------


  np=p_count;
  fvec(1:NDIM,1:np)=0.;
  Res(1:MAXCONN)=.false.
  
  do i=1,np
     
     do j=i+1,np
        lij=pos(:,i)-pos(:,j);
        l= sqrt( lij(1)*lij(1)  + lij(2)*lij(2)  + lij(3)*lij(3) );

        vij3=v(:,i)-v(:,j);
        vij=sqrt(vij3(1)*vij3(1) + vij3(2)*vij3(2) + vij3(3)*vij3(3));
        
        r_dot_v= (lij(1)*vij3(1) + lij(2)*vij3(2) + lij(3)*vij3(3))

        if ((l<=rc).and.(l > 0.)) then
           a_rep=Aij(btypes(i),btypes(j))
           
           ! Conservative 
           Fij_C= a_rep * (1.-l/rc)
           ! Dissipative
           Fij_D= 0.5 * sig_sq * (1.-l)*(1.-l)* (r_dot_v)  
           ! Random (brownian)
           call RAND_GEN(zeta)
           zeta=zeta-0.5;
           !write(*,*) 'zeta',zeta
           Fij_R= - sig * (1.-l) * zeta / sqrt_dt;
           ! Notice the minus sign when comparing to original paper (Groot and Rabone 2001)
           ! because the definition of rij is flipped                     
           fval= Fij_C + Fij_R + Fij_D ; 
           fvec(:,i)=fvec(:,i) + fval * lij/l;
           fvec(:,j)=fvec(:,j) - fval * lij/l;
        end if
        !write(*,*)parents(i),parents(j)
        !write(*,*)internalIndex(1:17)
        ! Bonded forces
        if ((parents(i)==parents(j)).and.(parents(i)<=pt_numBodies)) then
           BodyType = int(parentType(i));
           !write(*,*)BodyType
           bead1= int(internalIndex(i));
           bead2= int(internalIndex(j));
           !write(*,*) bead1, bead2,i,j
           !write(*,*)Connect(BodyType)%links(bead1,:)
           Res=(Connect(BodyType)%links(bead1,:).eq.bead2);
           !write(*,*)'Res=',Res
           !stop
           if (ANY(Res)) then
              fvec(:,i)=fvec(:,i) - 4.*lij;
              fvec(:,j)=fvec(:,j) + 4.*lij;   
           end if
        end if
     end do
     !stop
  end do
 
end subroutine pt_dpdNonBondedForces


!-------

subroutine rand_gen(zij)
  
  
  !***********************************
  !     RANDOM NUMBER GENERATOR
  !     HUSSEIN EZZELDIN UMCP FEB 2011
  !************************************
  implicit none
  
  INTEGER  ::   i_seed, dt_seed(8)
  INTEGER,ALLOCATABLE  :: a_seed(:)
  
  REAL(8),INTENT(OUT) :: zij
  
  i_seed = 8;
  
  CALL RANDOM_SEED(size=i_seed) 
  ALLOCATE(a_seed(1:i_seed)) 
  
  CALL RANDOM_SEED(get=a_seed) 
  
  CALL DATE_AND_TIME(values=dt_seed) 
  a_seed(i_seed)=dt_seed(1); 
  a_seed(1)=dt_seed(4)*dt_seed(7)*dt_seed(6) 
  
  CALL RANDOM_SEED(put=a_seed) 
  DEALLOCATE(a_seed) 
  
  CALL random_number(zij)
  
  
end subroutine rand_gen
