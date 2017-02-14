!!****if* source/Simulation/SimulationMain/WindTunnel/Simulation_defineDomain
!!
!! NAME
!!
!!  Simulation_defineDomain
!!
!! SYNOPSIS
!!
!!  Simulation_defineDomain(integer(MDIM), intent(in) :: nblks,
!!                          logical(:,:,:),intent(out) :: initialDomain,
!!                          integer(:,:,:,:)intent(out) ::boundaries)
!!                          
!!
!! DESCRIPTION
!!
!!  This routine provides a hook for defining domains that are not 
!!  shaped like boxes. The stub implementation returns a domain that is
!!  shaped like a  box. The users wishing to define a non-box like domain
!!  should have a custom implementation of this routine in their 
!!  Simulations directory.  Only domains that can be represented as
!!  a union of non-overlapping same-size boxes place on an equally-space
!!  grid can represented.
!!
!!  This feature is currently only available with PARAMESH. When building
!!  FLASH with UG, any Simulation_defineDomain is ignored.
!!
!!  Conceptually, one defines an irregular domain by placing obstacle blocks
!!  in a regular box-shaped domain. These obstacle blocks all have the
!!  same physical size, and the original domain size must be a multiple
!!  of this size in each dimesion. Obstacle may touch the outer boundaries
!!  and each other, and can be combined to form walls etc.
!!
!!  To use this feature:
!!  1) Set up a regular domain as usual. Use the nblockx, nblocky, and
!!     nblockz runtime parameters to divide the domain equally into
!!     a grid of  nblockx*nblocky*nblockz  original blocks.
!!  2) Provide an implementation of Simulation_defineDomain that marks
!!     some of these original blocks with obstacles. Only the parts of
!!     the original domain that are not so marked will end up being in
!!     the physical domain.
!!  PARAMESH will use the blocks that remain in the domain as root blocks
!!  (with refinement level 1) and may adaptively refine (but not combine)
!!  them as usual.
!!
!!  To implement Simulation_defineDomain: For each block that belongs in the domain
!!  the corresponding entry in the output array initialDomain is set true, and
!!  the corresponding entries for all faces in the boundaries array are set to
!!  "NOT_BOUNDARY". For the blocks that do not belong in the domain, the 
!!  corresponding initialDomain entries are set false. The boundaries entries 
!!  on faces that are adjacent to physical boundaries or to other blocks 
!!  not in the domain are redundant. However, those boundaries entries that 
!!  correspond to faces adjacent to the valid blocks in the domain must be
!!  given appropriate boundary values. The example included illustrates how
!!  to define a domain. 
!!
!!  This routine will be called with nblks set to the values given by the
!!  paramesh runtime parameters nblockx, nblocky, and nblockz, and the
!!  initialDomain and boundaries arrays will be dimensioned accordingly.
!!  That is, the boundaries set by this routine apply to blocks at refinement
!!  level 1. Any additional refinement levels are generated only after this
!!  routine returns.
!!
!! ARGUMENTS
!!
!!   nblks - Integer Array of size MDIM, containing the number of blocks along
!!           each dimension
!!   initialDomain - Logical 3D array of size nblks(IAXIS),nblks(JAXIS),
!!                                            nblks(KAXIS)
!!                   If the block i,j,k is in the domain then on return
!!                   initialDomain(i,j,k)=.true. otherwise .false.
!!   boundaries  - integer array of size 2*MDIM,nblks(IAXIS),nblks(JAXIS),
!!                                       nblks(KAXIS)
!!                   If the block i,j,k is in the domain, the values in
!!                   this array have no meaning. If the block is not in the
!!                   domain, but is adjacent to one or more blocks that
!!                   are in the domain, then the entries corresponding to
!!                   the common faces get initialized with appropriate
!!                   boundary values.
!!                   Constants for boundary conditions are defined in
!!                   "constants.h". Any boundary conditions in the allowed
!!                   range of numerical values -50..-20 can be used if they
!!                   are recognized by the GridBoundaryconditions implementation,
!!                   with the exception that PERIODIC is not allowed here.
!!
!! Example  
!!    Consider a 2D domain with the following shape with reflecting boundaries
!!        
!!        **** **** ****
!!       *              *
!!       *              *
!!       *               ****
!!       *                   *
!!       *                   *
!!       *                   *
!!       *                   *  
!!       *                   *  
!!       *                   *  
!!        ****               *  
!!            *              *        
!!            *              *        
!!             **** **** ****
!!
!!     It can be divided into 4x4 blocks as follows
!!        **** **** **** ****
!!       * T  * T  * T  * F  *
!!       *    *    *    *    *
!!        **** **** **** ****
!!       * T  * T  * T  * T  *
!!       *    *    *    *    *
!!        **** **** **** ****
!!       * T  * T  * T  * T  *
!!       *    *    *    *    *
!!        **** **** **** ****
!!       * F  * T  * T  * T  *
!!       *    *    *    *    *
!!        **** **** **** ****
!!
!! This domain can be initialize through the following sequence of
!!  instructions
!! 
!!   !! This implementation is meant to be used with NDIM==2.
!!   initialDomain = .true.
!!   boundaries = NOT_BOUNDARY
!!   initialDomain(1,1,1)=.false.
!!   boundaries(IHI_FACE,1,1,1)=REFLECTING
!!   boundaries(JHI_FACE,1,1,1)=REFLECTING
!!   initialDomain(4,4)=.false.
!!   boundaries(ILO_FACE,4,4,1)=REFLECTING
!!   boundaries(JLO_FACE,4,4,1)=REFLECTING
!! 
!! NOTES
!!
!!  To use any Simulation_defineDomain implementation that actually
!!  removes blocks with Paramesh 4, it is currently necessary to increase
!!  the value of the constant NBOUNDARIES. This preprocessor symbol is used
!!  by Paramesh 4 to dimension certain arrays that store information on
!!  boundary blocks. The value of NBOUNDARIES should be set to
!!    2*NDIM  +  (number of blocks for which this routine returns false).
!!  This is most conveniently done with a PPDEFINE in a Config file before
!!  setup. See the SodStep simulation as an example. FLASH will abort if
!!  Simulation_defineDomain is used and NBOUNDARIES is too low.
!!
!!  For users' convenience, the constants defining the low 
!!  and high faces along the dimensions are defined in "constants.h".
!!  They can be used to specify the first index of the boundaries
!!  array. They are {I,J,K}LO_FACE and {I,J,K}HI_FACE, respectively.
!!
!!  Constants for boundary conditions like REFLECTING and OUTFLOW
!!  are defined in "constants.h" as well. The most common of them
!!  correspond to recognized values of the {xyz}{lr}_boundary_type
!!  runtime parameters.
!!
!!  This routine is only meaningful (and will only be called) when an
!!  AMR Grid implementation is used.
!!
!!  Particles have not been tested in combination with arbitrary
!!  boundaries. Particles (active or passive) should not be used
!!  in irregular domain, at least not if there is any chance that
!!  a particle will run into an obstacle block. The behavior of
!!  particles at obstacle boundaries should be cinsidered as
!!  undefined.
!!
!!***

subroutine Simulation_defineDomain(initialDomain,boundaries,nblks)
#include "constants.h"
  implicit none
  integer,dimension(MDIM),intent(IN) :: nblks
  integer,dimension(2*MDIM,nblks(IAXIS),nblks(JAXIS),nblks(KAXIS)),&
       intent(OUT)::boundaries
  logical,dimension(nblks(IAXIS),nblks(JAXIS),nblks(KAXIS)),&
       intent(OUT)::initialDomain

  integer :: x_sweep

  initialDomain =.true.
  boundaries=NOT_BOUNDARY

! Leave 3 blocks just in front of the step.
  do x_sweep = 4, nblks(IAXIS)
     initialDomain(x_sweep, 1, 1) = .false.
     boundaries(:,x_sweep,1,1)=REFLECTING
  end do


end subroutine Simulation_defineDomain
