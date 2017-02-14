!!****if* source/diagnostics/ProtonImaging/ProtonImagingMain/ProtonBeams/pi_capsuleNextGrainIndices
!!
!! NAME
!!
!!  pi_capsuleNextGrainIndices
!!
!! SYNOPSIS
!!
!!  call pi_capsuleNextGrainIndices (integer, intent (in)    :: L,
!!                                   integer, intent (inout) :: i,
!!                                   integer, intent (inout) :: j,
!!                                   integer, intent (inout) :: k,
!!                                   logical, intent (out)   :: valid)
!!
!! DESCRIPTION
!!
!!  This routine retrieves the next proper capsule grain indices from a (presumably) valid
!!  capsule grain index triple (i,j,k). Even if the input (i,j,k) does not correspond to
!!  a valid capsule grain index triple, the routine finds the next valid capsule grain index
!!  triple.
!!
!!  The grain level L is synonymous to the number of cubes with side lengths (R/L) placed
!!  along the radius R of the capsule sphere. Each cube constitutes one particular grain
!!  of the capsule and is given a labelling index triple (the capsule grain indices).
!!  The following problem setup explains how these cubes are labeled and the diophantine
!!  equation that defins them. 
!!
!!  Problem setup:
!!  -------------
!!
!!  The following setup of the problem is assumed:
!!
!!      1) A 3D cartesian coordinate system is placed in space.
!!      1) A sphere is placed at the cartesian origin (0,0,0), having exactly L equal sized
!!         divisions along its radius R in all three cartesian directions.
!!      2) Cubes of side size R/L are placed inside the sphere with one vertex of each of
!!         the 8 innermost cubes located at the origin (0,0,0) of the sphere and the cartesian
!!         coordinate system.
!!      3) Each cube (capsule grain) can be labeled by three indices (i,j,k), where (i,j,k)
!!         denotes the largest number combination of the vertex indices. A vertex index simply
!!         denotes the coordinate number on which it is placed in the coordiate system.
!!
!!  As an example, consider the innermost cube with its 8 vertex indices (0,0,0), (0,0,1),
!!  (0,1,0), (1,0,0), (0,1,1), (1,0,1), (1,1,0) and (1,1,1). The cube's labelling indices would
!!  then be the last vertex indices (1,1,1).
!!
!!  Each cube's labelling indices are therefore an indication of how far out it reaches in space.
!!  To determine if a cube is located inside the sphere, its labelling indices (i,j,k) must
!!  obey the following inequality:
!!
!!                        (i*R/L)^2 + (j*R/L)^2 + (k*R/L)^2 <= R^2
!!
!!  which (after rescaling to unit radius of the sphere) is the same as:
!!
!!                              i^2 + j^2 + k^2 <= L^2       (I)
!!
!!  The present routine finds, from a (valid or not) index triple the next index triple,
!!  such that the diophantine equation (I) is satisfied.
!!
!!  The algorithm used is quite simple. It is an open loop that checks the i,j,k indices
!!  between -L and + L and excludes those that violate equation (I) and all those (i,j,k)
!!  triples that contain at least a zero index. Index i varies fastest, index k slowest.
!!
!!
!! ARGUMENTS
!!
!!  L     : the capsule grain level
!!  i     : i-index of the original (valid or not) and new triple (i,j,k).
!!  j     : j-index of the original (valid or not) and new triple (i,j,k).
!!  k     : k-index of the original (valid or not) and new triple (i,j,k).
!!  valid : if true, a valid new index triple has been found. False, if no more triples possible. 
!!
!!***

subroutine pi_capsuleNextGrainIndices (L,i,j,k,valid)

  implicit none

  integer, intent (in)    :: L
  integer, intent (inout) :: i,j,k
  logical, intent (out)   :: valid
!
!
!     ...Find the next index triple.
!
!
  i = i + 1

  do while (     i * i + j * j + k * k > L * L &    ! check violation of equation (I)
            .or. i > L - 1                     &    ! largest value for an index is L - 1
            .or. i == 0                        &    ! 
            .or. j == 0                        &    ! exclude zeros
            .or. k == 0                        )    !

     i = i + 1
     if (i > L - 1) then
         i = - L + 1
         j = j + 1
         if (j > L - 1) then
             j = - L + 1
             k = k + 1
             if (k > L - 1) then
                 valid = .false.
                 return
             end if
         end if
     end if
  end do

  valid = .true.
!
!
!     ...Ready!
!
!
  return
end subroutine pi_capsuleNextGrainIndices
