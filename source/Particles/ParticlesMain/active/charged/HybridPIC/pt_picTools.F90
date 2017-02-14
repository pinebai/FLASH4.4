!!****ih* source/Particles/ParticlesMain/active/charged/HybridPIC/pt_picTools
!!
!! NAME
!!  pt_picTools
!!
!! SYNOPSIS
!!
!!  use pt_picTools
!!
!! DESCRIPTION 
!!  
!!  Module containing general particle routines.
!!  Independent of other FLASH routines and data.
!!  Does not call anything problem specific or FLASH related,
!!  except for ParticleAdd and ParticleDelete.
!!***

module pt_picTools
  implicit none

#include "Flash.h"
#include "constants.h"  
#include "Particles.h"

  integer, save :: n_proc       ! number of processors

contains

  subroutine pt_picBoundaryVelocity(v, n, u, t, m)
    ! Generate a random particle velocity V [m/s], for a particle emerging 
    ! from a boundary with normal N.  The boundary is moving with velocity 
    ! U [m/s].
    ! The particle distribution has temperature t [K], and the particles has 
    ! mass M [kg]. 
    ! Note that the distribution is skewed, following Garcia p. 369. 

    implicit none
    real, intent(out) :: v(3)
    real, intent(in)  :: n(3), u(3)
    real, intent(in)  :: t, m

    ! Should get this from FLASH
    real, parameter   :: pi = 3.141592653589793238462643383279503
    real, parameter   :: kb = 1.38066e-23 ! Boltzmann [J/K = kg(m/s)^2/K]
    real :: nx(3), ny(3), nz(3)
    real :: a, pw, vmax, maxp, pr, w(3)

    a = m/(2*kb*t)
    pw = sqrt(1/a)         ! approx. width of distribution

    call pt_picRan_exp(w(1))
    w(1) = sqrt(w(1))*pw
    call pt_picRan_normal(w(2))
    w(2) = w(2)*pw/sqrt(2.)
    call pt_picRan_normal(w(3))
    w(3) = w(3)*pw/sqrt(2.)

    ! W is now a random velocity with its x-component along N
    nx = n/sqrt(sum(n*n))         ! normalized normal
    ! find a vector orthogonal to nx
    if (abs(nx(1)) .lt. .9) then
       ny(1) = 1; ny(2) = 0; ny(3) = 0
    else
       ny(1) = 0; ny(2) = 1; ny(3) = 0
    end if
    ny = ny - sum(ny*nx)*nx
    ny = ny/sqrt(sum(ny*ny))
    ! nz = nx x ny
    nz(1) = nx(2)*ny(3)-nx(3)*ny(2)
    nz(2) = nx(3)*ny(1)-nx(1)*ny(3)
    nz(3) = nx(1)*ny(2)-nx(2)*ny(1)
    ! now (nx,ny,nz) is a basis for w
    v = w(1)*nx + w(2)*ny + w(3)*nz + u   ! v in the simulation frame
  end subroutine pt_picBoundaryVelocity

  subroutine pt_picMaxwellian(u, v, t, m)
    ! U is a random velocity from the Maxwellian with bulk velocity V, 
    ! and temperature T for a species of mass M. 
    implicit none
    real, intent(out) :: u(3)         ! [m/s]
    real, intent(in)  :: v(3), t, m   ! [m/s], [K], [kg]

    real :: a, pw
    real, parameter   :: kb = 1.38066e-23 ! Boltzmann [J/K = kg(m/s)^2/K]

    a = m/(2*kb*t)
    pw = sqrt(1/a)         ! approx. width of distribution

    call pt_picRan_normal(u(1))
    call pt_picRan_normal(u(2))
    call pt_picRan_normal(u(3))
    u = u*pw/sqrt(2.)

    u = u + v
  end subroutine pt_picMaxwellian

  subroutine pt_picSpheran(r)
    ! Based on undervisning/prog/tentor/990830/spheran.c
    ! Computes the rectangular coordinates of a random point on the 
    ! unit sphere.  Algoritm from Knuth 

    implicit none
    real, intent(out) :: r(3)
    real :: u(2), a, s

    do 
       call pt_picRan_uni(u(1))
       call pt_picRan_uni(u(2))
       u = u*2-1                  ! [-1,1]
       s = sum(u*u)
       if (s <= 1) exit
    end do
    a = 2*sqrt(1-s)

    r(1:2) = a*u
    r(3) = 2*s-1
  end subroutine pt_picSpheran

  real function angle(u, v)
    ! Should be general... see simf.h
    implicit none
    real, intent(in)  :: u(3)
    real, intent(in)  :: v(3)
    real :: denom

    denom = sqrt(sum(u*u)*sum(v*v))
    if (denom == 0.0) print *, 'Error.  Division by zero in angle()'

    angle = acos(sum(u*v)/denom)
  end function angle

  subroutine pt_picRec2sph(d, lng, lat, r)
    ! Convert rectangular coordinates to spherical.
    ! r=(1,0,0) => 1, 0, 0;  r=(0,1,0) => 1, pi/2, 0;  r=(0,0,1) => 1, 0, pi/2
    ! longitude in [-pi,pi] and latitude in [-pi/2,pi/2]
    ! See also mercury/mc/src/mercion/merctraj.cpp
    implicit none
    real, intent(out) :: d, lng, lat
    real, intent(in)  :: r(3)
    real              :: s(3)
    real, parameter   :: pi = 3.1415926535897932

    d = sqrt(sum(r*r))
    s = r/d                     ! normalized r
    lat = acos(s(3))            ! [0,pi] for z=1 to z=-1
    lat = pi/2 - lat            ! [-pi/2,pi/2]
    lng = atan2(s(1), s(2))     ! (-pi,pi]
  end subroutine pt_picRec2sph

  subroutine pt_picRan_ini(rpSeed, n, np)
    use ut_randomInterface, only: ut_randomSeed, ut_randomNumber
    implicit none
    integer, intent(in) :: rpSeed, n, np
    integer :: i, seedSize
    real :: u
    integer,allocatable :: seedArray(:)

    if (rpSeed .GE. 0) then
       call ut_randomSeed(ut_size=seedSize)
       allocate(seedArray(seedSize))
       seedArray(1:seedSize) = rpSeed
       call ut_randomSeed(ut_put=seedArray)
       deallocate(seedArray)
    end if
    n_proc = np
    do i = 1, n              ! put all processors at different numbers
       call ut_randomNumber(u)
    end do
  end subroutine pt_picRan_ini


  subroutine pt_picRan_uni(u)
    ! http://en.wikipedia.org/wiki/Mersenne_twister
    use ut_randomInterface, only: ut_randomNumber
    implicit none
    real, intent(out) :: u
    integer :: i

    do i = 1, n_proc            ! all processors get different numbers
       call ut_randomNumber(u)
    end do
  end subroutine pt_picRan_uni

  integer function pt_picRan_int(l, u)
    ! Random integer in [l,u]
    implicit none
    integer, intent(in) :: l, u
    real :: r
    call pt_picRan_uni(r)
    pt_picRan_int = floor(l + r*(u-l) + 0.5) ! OK???
  end function pt_picRan_int

  subroutine pt_picRan_exp(u)
    implicit none
    real, intent(out) :: u
    real              :: r
    do
       call pt_picRan_uni(r)
       if (r .gt. 0.0) exit
    end do
    u = -log(r)
  end subroutine pt_picRan_exp

  subroutine pt_picRan_expr(u, r)
    implicit none
    real, intent(out) :: u
    real, intent(in)  :: r   ! Rate [s^{-1}]
    real              :: x

    call pt_picRan_exp(x)
    u = x/r
  end subroutine pt_picRan_expr

  subroutine pt_picRan_normal(u)
    ! Gaussian distribution by Box-Muller. 
    ! We are wasting one random nmber, see 001219/boxmuller.c
    ! Also, it can be done without trig.  See NR. 
    implicit none
    real, intent(out) :: u
    real, parameter   :: pi = 3.1415926535897932
    real              :: r1, r2
    call pt_picRan_uni(r1)
    do
       call pt_picRan_uni(r2)
       if (r2 .gt. 0.0) exit
    end do
    u = sqrt(-2*log(r2))*sin(2*pi*r1)
  end subroutine pt_picRan_normal


!!****if* hybrid/pt_picTools.F90
!! NAME
!!   pt_picParticleAdd -- Adds a new particle given its properties
!! SYNOPSIS
!!   call pt_picParticleAdd(r, v, mass, charge, specie, blknr)
!! FUNCTION
!!   Adds a new particle for the local processor with specified properties. 
!!   Requires that the position, r, is in the local block, blknr. 
!!   Also requires that the block is a leaf. 
!!   After adding particles, ModuleReStuffArray::ReStuffArray needs 
!!   to be called before iterating over the particles again so that 
!!   the particle array is compacted and tree::lnblocks is updated. 
!! INPUTS
!!   r        - The 3D position of the particle 
!!   v        - The 3D velocity of the particle 
!!   mass     - The mass of the particle
!!   charge   - The charge of the particle
!!   specie   - An integer number denoting the species of the particle, e.g., 
!!              1 for protons. 
!!   blknr    - The integer number of the local block that the particle 
!!              position, r, falls inside.  1 <= blknr <= tree::lnblocks
!! NOTES
!!   We could do without the block parameter and instead call 
!!   ClassifyParticleLocation, FindParticleBlock in ParticleMeshInfo.F90 
!!   to get the block number.  This does not work however.  Also, this 
!!   only makes sense if all processes produces all particles. 
!! AUTHOR
!!  Mats Holmstrom, mailto:matsh@irf.se
!! CREATION DATE
!!  January 2007
!!***
  subroutine pt_picParticleAdd(r, v, mass, charge, specie, blknr)
    use Particles_data, only: particles, pt_numLocal, pt_maxPerProc
    
    implicit none
    real, intent(in)    :: r(3), v(3), mass, charge
    integer, intent(in) :: specie, blknr
    logical, save :: first_call = .true.
    integer :: i

    if (first_call) then 
       particles(BLK_PART_PROP,:) = NONEXISTENT
       first_call = .false.
    end if

    if (pt_numLocal >= pt_maxPerProc) then
       print *, 'pt_initPositions:  max part/proc = ', pt_maxPerProc
    else
       pt_numLocal = pt_numLocal + 1
       i = pt_numLocal
       ! insert the new particle at the end
       particles(BLK_PART_PROP, i) = blknr
       particles(TAG_PART_PROP, i) = 0       
! See ./Particles/ParticlesInitialization/pt_createTag.F90 for tag generation
! *** No way to generate one tag?  Need to generate for ALL particles? ***

       particles(MASS_PART_PROP, i) = mass
       particles(CHARGE_PART_PROP, i) = charge
       particles(SPECIE_PART_PROP, i) = specie
       ! Particle position and velocity.
       particles(POSX_PART_PROP, i) = r(1)
       particles(VELX_PART_PROP, i) = v(1)
       particles(POSY_PART_PROP, i) = r(2)
       particles(VELY_PART_PROP, i) = v(2)
       particles(POSZ_PART_PROP, i) = r(3)
       particles(VELZ_PART_PROP, i) = v(3)
    endif
end subroutine pt_picParticleAdd


!!****if* hybrid/pt_picTools.F90
!! NAME
!!   pt_picParticleDelete -- Deletes a specified particle 
!! SYNOPSIS
!!   call pt_picParticleDelete(i)
!! FUNCTION
!!   Requires that the particle exist. 
!!   After deleting particles, ModuleReStuffArray::ReStuffArray needs 
!!   to be called before iterating over the particles again so that 
!!   the particle array is compacted and tree::lnblocks is updated. 
!! INPUTS
!!   i   - The integer index of the particle that will be deleted 
!! NOTES
!!   The deletion is done as in ParticleBoundaries::Outflow
!!   CANNOT be uded in loop bounded by pt_numLocal since we modify it!
!! AUTHOR
!!  Mats Holmstrom, mailto:matsh@irf.se
!! CREATION DATE
!!  January 2007
!!***
  subroutine pt_picParticleDelete(i)
    use Particles_data, only: particles, pt_numLocal
    implicit none
    integer, intent(in) :: i

    if (particles(BLK_PART_PROP, i) == NONEXISTENT) then
       print *, 'pt_picParticleDelete(', i, ') NONEXISTENT. Impossible...'
    else
       ! Replace with last particle, as in gr_ptHandleExcess.F90
       particles(1:NPART_PROPS, i) = particles(1:NPART_PROPS, pt_numLocal)
       ! Make the particle invalid, as in gr_ptMoveOffProc.F90
       particles(BLK_PART_PROP, pt_numLocal) = NONEXISTENT
       pt_numLocal = pt_numLocal-1
    end if
  end subroutine pt_picParticleDelete
end module pt_picTools
