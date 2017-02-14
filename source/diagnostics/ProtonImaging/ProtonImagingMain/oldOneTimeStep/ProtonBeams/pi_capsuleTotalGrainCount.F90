!!****if* source/diagnostics/ProtonImaging/ProtonImagingMain/oldOneTimeStep/ProtonBeams/pi_capsuleTotalGrainCount
!!
!! NAME
!!
!!  pi_capsuleTotalGrainCount
!!
!! SYNOPSIS
!!
!!  pi_capsuleTotalGrainCount (integer, intent (in) :: L)
!!
!! DESCRIPTION
!!
!!  This function gives the total number of capsule grains for a specific capsule grain level L.
!!  The grain level L is synonymous to the number of cubes with side lengths (R/L) placed
!!  along the radius R of the capsule sphere. The total grain count inside the capsule is
!!  thus equal to the total number of cubes with side lengths (R/L) that can be completely
!!  placed inside the spherical capsule, if the center of the capsule coincides with a vertex
!!  of the innermost 8 cubes. The following problem setup explains how these cubes are labeled
!!  and the diophantine equation that needs to be solved. 
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
!!      3) Each cube can be labeled by three indices (i,j,k), where (i,j,k) denotes the
!!         largest number combination of the vertex indices. A vertex index simply denotes the
!!         coordinate number on which it is placed in the coordiate system.
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
!!  The present function counts the total number of integer solutions of the diophantine
!!  equation (I), where |i,j,k| >= 1. Note that zeros for the i,j,k are excluded.
!!
!!  Since for L = 1 there would be no grains in the capsule but we still want to have one
!!  grain located at the center of the capsule, the L = 1 case is treated separately and
!!  given a count of 1. The index triple returned in this case is (0,0,0).
!!
!!  The following table gives an indication on how many grains are found for each grain
!!  level L:
!!
!!                              L | number of grains
!!                             ---------------------
!!                              1 |           1
!!                              2 |           8
!!                              3 |          56
!!                              4 |         136
!!                              5 |         304
!!                             10 |        3280
!!                             20 |       29752
!!                            100 |     4094208
!!                            200 |    33132200
!!
!!  Algorithmic considerations:
!!  --------------------------
!!
!!  Since for large L the result is related to Pi = 3.14159...., there is no simple algebraic
!!  expression for this number. Rather than applying a triple 'brute force' loop between -L
!!  and +L for each index and checking the validity of equation (I), we can do better by
!!  determining largest cubes and squares inside specific volumes and areas. Index count in
!!  cubes and squares is simple and straightforward. The remaining areas near the sphere's
!!  surface area are the only ones that need explicit checking of equation (I). Oh-group
!!  symmetry of the capsule sphere with its inserted cubes is exploited and hence the check
!!  for valid (i,j,k) indices is reduced to one octant only.
!!
!! ARGUMENTS
!!
!!  L : the capsule grain level
!!
!! NOTES
!!
!!  Returns a value of -1, if the total number of grains is expected to be larger than
!!  the largest integer representable on the current machine.
!!
!!***

integer function pi_capsuleTotalGrainCount (L)

  implicit none

  integer, intent (in) :: L

  integer :: f,i,j,k,m,p,q,r
  integer :: kCube, jSquare

  real    :: invSqrt3
!
!
!     ...Handle special case of grain level = 1.
!
!
  if (L == 1) then
      pi_capsuleTotalGrainCount = 1
      return
  end if
!
!
!     ...As a precaution, we need to estimate, if the number of total grains will be
!        larger than the largest integer representable on the current machine. If this
!        is the case, the function returns a value of -1, which serves as a handle to
!        indicate the overflow problem.
!
!        For large L, the approximate number of cubes in the sphere approaches the
!        volume formula for the sphere: (4Pi/3) * L^3 from above (i.e. the volume
!        formula gives larger value than the true number of cubes). Hence, overflow
!        of the volume formula can be used as a limiting factor for L. We have:
!
!                    (4Pi/3) * L^3  >  LPI     (LPI = Largest Positive Integer)
!
!        for the condition of L. This can be transformed to:
!
!                           L > ((3/4Pi) * LPI) ^ (1/3)
!
!        Note, that if (3/4Pi) is only approximated from below, then the critical L
!        will be slightly lower and we are safe. Here we take (3/4Pi) = 0.238732 as
!        a lower bound for the true value 0.2387324146....
!
!
  if (L > int ((0.238732 * huge (1)) ** (1.0/3.0))) then
      pi_capsuleTotalGrainCount = -1
      return
  end if
!
!
!     ...We are good to go. First, determine largest enclosed cube and
!        count all cubes inside this large cube.
!
!
  invSqrt3 = 1.0 / sqrt (3.0)

  m = 0
  p = L * L

  kCube = int (L * invSqrt3) + 1  ! get oversized cube
  q = kCube * kCube
  do while (q + q + q >= p)        ! find proper inside cube
     kCube = kCube - 1
     q = kCube * kCube
  end do

  m = m + 8 * q * kCube           ! all indices in cube of side length 2 * kCube
  p = p - q
!
!
!     ...Next, for one octant cube of the large cube, determine the largest
!        enclosed square (of cubes). There are 3 such identical squares on the
!        outer sides of the octant cube and there are 8 octant cubes altogether
!        => a factor of 24 is needed for each square.
!
!
  do k = kCube + 1, L - 1
     p = p - k - k + 1

     jSquare = int (sqrt (real (p) * 0.5)) + 1    ! get oversized square
     q = jSquare * jSquare
     do while (q + q >= p)                        ! find proper inside square
        jSquare = jSquare - 1
        q = jSquare * jSquare
     end do

     m = m + 24 * q                               ! each square is added to 3 cube sides per octant
     q = p - q
!
!
!     ...The remaining cubes, flanking both sides of the current square. An additional
!        factor of 2 is needed, except for those cubes which have equal indices (j=k).
!
!
     do j = jSquare + 1, k
        q = q - j - j + 1
        r = q
        f = 1 + min (1,k-j)
        do i = 1, j - 1
           r = r - i - i + 1
           if (r >= 0) then
               m = m + 24 * f
           else
               exit
           end if
        end do
     end do

  end do     ! next square

  pi_capsuleTotalGrainCount = m
!
!
!     ...Ready!
!
!
  return
end function pi_capsuleTotalGrainCount
