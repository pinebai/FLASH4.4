!!***h* source/physics/ImBound/ImBoundMain/ImBound.h
!!
!! This is the internal header file for the 
!! Immersed Boundary Method Module
!!
!!***

#define IB_OPTIMIZE 1
!#define CHK_BND_VLC 1

!#define TEST_IB_FORCES 1
!!#define INVERSE_GCELL_FILL 1

!! Compute the tangent forces from vorticity field.
!!#define TANGENT_WITH_VORTICITY 1

!! In case tangent forces not computed from vorticity field,
!! if USE_CF defined, the tangent forces will be computed from 
!! the velocity difference along the streamline direction next
!! to the surface.
#define USE_CF 1

!! Case of computing tangent distributed forces from strain rate tensor.
#ifndef TANGENT_WITH_VORTICITY
#define IB_GET_DERIVS 1
#endif


#define FORCE_FLOW 1
#define COMPUTE_FORCES 2


#define TWO_NODE_SEGMENT_VERT  21
#define TWO_NODE_SEGMENT_CEN   22
#define THREE_NODE_TRIANG_VERT 31
#define THREE_NODE_TRIANG_CEN  32
#define FOUR_NODE_QUAD_VERT    41
#define FOUR_NODE_QUAD_CEN     42


