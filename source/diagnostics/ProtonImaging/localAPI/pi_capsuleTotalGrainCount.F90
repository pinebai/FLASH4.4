!!****if* source/diagnostics/ProtonImaging/localAPI/pi_capsuleTotalGrainCount
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

  pi_capsuleTotalGrainCount = 0

  return
end function pi_capsuleTotalGrainCount
