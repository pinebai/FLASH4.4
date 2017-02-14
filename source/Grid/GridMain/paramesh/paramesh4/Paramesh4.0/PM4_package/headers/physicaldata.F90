!----------------------------------------------------------------------
! PARAMESH - an adaptive mesh library.
! Copyright (C) 2003
!
! Use of the PARAMESH software is governed by the terms of the
! usage agreement which can be found in the file
! 'PARAMESH_USERS_AGREEMENT' in the main paramesh directory.
!----------------------------------------------------------------------

!!****h* headers/physicaldata
!!
!! NAME
!!
!!   physicaldata
!! 
!! SYNOPSIS
!!
!!   module physicaldata
!!  
!! USES
!!
!!   paramesh_dimensions
!!
!! DESCRIPTION
!!
!!   Fortran 90 Module which 'holds' the main solution data structures of
!!   PARAMESH.  These include 'unk', 'facevarx', 'facevary', 'facevarz', 
!!   'unk_e_x', 'unk_e_y', 'unk_e_z', and 'unk_n'.  This module also holds the
!!   data used in by the flux conservation routines (i.e. flux_x, flux_y, and 
!!   flux_z).  Other data which controls how the solution data arrays are 
!!   manipulated by PARAMESH are also kept in this module such as the control 
!!   arrays for selecting the algorithms used for interpolation on guardcells 
!!   filling (i.e. interp_mask_*) and for selective guardcell filling (i.e. 
!!   gcell_on_*).  Arrays storing geometry information (i.e. cell volumes, 
!!   cell-face areas etc.) are also kept here.
!!
!! AUTHORS
!!
!!  Peter MacNeice and Kevin Olson
!!
!!***


!!REORDER(5): unk, facevar[xyz], tfacevar[xyz]
!!REORDER(4): recvar[xyz]f
#include "paramesh_preprocessor.fh"


!------------------------------------------------------------------------------
! Physicaldata module
!------------------------------------------------------------------------------




      module physicaldata

      use paramesh_dimensions

      private


!----------------------
! Solution Variables
!----------------------

!---------------------------------------
! Allocate memory for solution variables
!---------------------------------------

!------------------------------------------------------------------------------
! CELL CENTERED CENTERED DATA
!------------------------------------------------------------------------------


! the solution for cell-centered quantities.
      public :: unk, interp_mask_unk, interp_mask_unk_res
      public :: gcell_on_cc,int_gcell_on_cc
      public :: ngcell_on_cc
      public :: checkp_on_cc
      public :: gcell_on_cc_pointer
#ifndef LIBRARY
      real,save ::  unk(nvar,                                           &
     &                  il_bnd:iu_bnd,jl_bnd:ju_bnd,kl_bnd:ku_bnd,      &
     &                  maxblocks)
      integer, save :: interp_mask_unk(nvar)
      integer, save :: interp_mask_unk_res(nvar)
      integer, save :: gcell_on_cc_pointer(nvar)
      logical, save :: gcell_on_cc(nvar)
      logical, save :: int_gcell_on_cc(nvar)
      logical, save :: checkp_on_cc(nvar)
#else
      real,allocatable,save ::  unk(:,:,:,:,:)
      integer,allocatable,save :: interp_mask_unk(:)
      integer,allocatable,save :: interp_mask_unk_res(:)
      integer,allocatable,save :: gcell_on_cc_pointer(:)
      logical,allocatable,save :: gcell_on_cc(:)
      logical,allocatable,save :: int_gcell_on_cc(:)
      logical,allocatable,save :: checkp_on_cc(:)
#endif
      integer, save :: ngcell_on_cc
      target :: unk


!!****v* physicaldata/unk
!!
!! NAME
!!
!!   unk
!! 
!! SYNOPSIS
!!
!!   public, real, target :: 
!!    unk(nvar,il_bnd:iu_bnd,jl_bnd:ju_bnd,kl_bnd:ku_bnd,maxblocks)
!!
!!   public, real, allocatable, target :: unk(:,:,:,:,:)
!!
!! DESCRIPTION
!!
!!   Array for holding cell centered data.  May or may not have permanent 
!!   guardcell storage allocated. If the LIBRARY flag is defined 
!!   paramesh_preprocessor.fh this array is allocated at runtime.
!!
!!***

!!****v* physicaldata/interp_mask_unk
!!
!! NAME
!!
!!   interp_mask_unk
!! 
!! SYNOPSIS
!!
!!   public, integer :: interp_mask_unk(nvar)
!!
!!   public, integer, allocatable :: interp_mask_unk(:)
!!
!! DESCRIPTION
!!
!!   Array to select different interpolation orders (or other user defined 
!!   interpolation schemes) used during prolongation operations on the cell 
!!   centered data stored in unk.  The values stored in interp_mask_unk may 
!!   be changed at anytime during execution and they can take on different 
!!   values for each different variable stored in unk.
!!
!!   The default value(s) for interp_mask_unk is 1 for all variables, 
!!   1 through nvar.
!!   The default behaviour for using different values of unk is:
!!     interp_mask_unk =  0  -> direct injection
!!     interp_mask_unk =  1  -> linear interpolation
!!     interp_mask_unk =  2  -> Lagrange polynomial of 2nd order
!!     interp_mask_unk =  3  -> Lagrange polynomial of 3rd order
!!     interp_mask_unk =  4  -> Lagrange polynomial of 4th order
!!     interp_mask_unk >= 20 -> user defined
!!
!!   Usage Example:
!!     If you wish to use the default interpolation during prolongation for 
!!     all variables in unk, but to use a 2nd order polynomial interpolation 
!!     for variable 4 then set,
!!
!!     interp_mask_unk(4) = 2
!!
!!***

!!****v* physicaldata/interp_mask_unk_res
!!
!! NAME
!!
!!   interp_mask_unk_res
!! 
!! SYNOPSIS
!!
!!   public, integer :: interp_mask_unk_res(nvar)
!!
!!   public, integer, allocatable :: interp_mask_unk_res(:)
!!
!! DESCRIPTION
!!
!!   Array to select different interpolation orders (or other user defined 
!!   interpolation schemes) used during restriction operations on the cell 
!!   centered data stored in unk.  The values stored in interp_mask_unk_res 
!!   may be changed at anytime during execution and they can take on different 
!!   values for each different variable stored in unk.
!!
!!   The default value(s) for interp_mask_unk_res is 1 for all variables, 
!!     1 through nvar.
!!   The default behaviour for using different values of unk is:
!!     interp_mask_unk_res =  0  -> illegal value
!!     interp_mask_unk_res =  1  -> Simple average of child data to parent 
!!                                  (linear interpolation)
!!     interp_mask_unk_res =  2  -> Lagrange polynomial of 2nd order
!!     interp_mask_unk_res =  3  -> Lagrange polynomial of 3rd order
!!     interp_mask_unk_res =  4  -> Lagrange polynomial of 4th order
!!     interp_mask_unk_res >= 20 -> user defined
!!   
!!   Usage Example:
!!     If you wish to use the default interpolation during restriction for 
!!     all variables in unk, but to use a 2nd order polynomial interpolation 
!!     for variable 4 then set,
!!
!!     interp_mask_unk_res(4) = 2
!!
!!***

!!****v* physicaldata/gcell_on_cc
!!
!! NAME
!!
!!   gcell_on_cc
!! 
!! SYNOPSIS
!!
!!   public, logical :: gcell_on_cc(nvar)
!!
!!   public, logical, allocatable :: gcell_on_cc(:)
!!
!! DESCRIPTION
!!
!!   A masking array used to control which variables in the unk array have 
!!   their guardcells filled during a call to amr_guardcell or 
!!   amr_1blk_guardcell.  If gcell_on_cc(ivar) is set to .true. then the 
!!   guardcells for variable 'ivar' in unk will be filled, if it is set to 
!!   .false. then variable 'ivar' will be skipped and its guardcells will not 
!!   be filled.  
!!
!!   The values in gcell_on_cc can be changed at any point during a run and 
!!   allow a user to avoid unecessary data movement in cases when only a small
!!   set of the variables in unk need guardcell values at some point during 
!!   their algorithm.
!!
!!***

!!****iv physicaldata/int_gcell_on_cc
!!
!! NAME
!!
!!   int_gcell_on_cc
!! 
!! SYNOPSIS
!!
!!   public, logical :: int_gcell_on_cc(nvar)
!!
!!   public, logical, allocatable :: int_gcell_on_cc(:)
!!
!! DESCRIPTION
!!
!!   ??? I don't know ???
!!
!!***

!!****iv physicaldata/ngcell_on_cc
!!
!! NAME
!!
!!   ngcell_on_cc
!! 
!! SYNOPSIS
!!
!!   public, integer :: ngcell_on_cc
!!
!! DESCRIPTION
!!
!!   An integer variable which holds the current number of variables in unk
!!   selected for guardcell filling by the mask array 'gcell_on_cc'.
!!
!!***

!!****iv physicaldata/gcell_on_cc_pointer
!!
!! NAME
!!
!!   gcell_on_cc_pointer
!! 
!! SYNOPSIS
!!
!!   public, integer :: gcell_on_cc_pointer(nvar)
!!
!!   public, integer, allocatable :: gcell_on_cc_pointer(:)
!!
!! DESCRIPTION
!!
!!   I don't know ????
!!
!!***

!!***v* physicaldata/checkp_on_cc
!!
!! NAME
!!
!!   checkp_on_cc
!! 
!! SYNOPSIS
!!
!!   public, logical :: checkp_on_cc(nvar)
!!
!!   public, logical, allocatable :: checkp_on_cc(:)
!!
!! DESCRIPTION
!!
!!   A logical mask array which controls which variables are added to a 
!!   checkpoint by calling the 'amr_checkpoint_wr' routine.  If 
!!   checkp_on_cc(ivar) is set to to true, then the ivar'th variable in the 
!!   unk array is added to the checkpoint file.  The default behaviour has ALL
!!   the values in checkp_on_cc set to .true. so that all variables in unk are
!!   written to the checkpoint file.
!!
!!***

!------------------------------------------------------------------------------
! FACE CENTERED DATA
!------------------------------------------------------------------------------


! the solution for cell-face-centered quantities.
      public :: facevarx,facevary,facevarz
      public :: interp_mask_facex,interp_mask_facey,interp_mask_facez
      public :: interp_mask_facex_res,interp_mask_facey_res, & 
     &          interp_mask_facez_res
      public :: gcell_on_fc,int_gcell_on_fc
      public :: ngcell_on_fc
      public :: gcell_on_fc_pointer
      public :: checkp_on_fc
#ifndef LIBRARY
      real,save ::  facevarx(nbndvar,                                   &
     &                       il_bnd:iu_bnd+1,jl_bnd:ju_bnd,             &
     &                       kl_bnd:ku_bnd,                             &
     &                       maxblocksf)
      real,save ::  facevary(nbndvar,                                   &
     &                       il_bnd:iu_bnd,jl_bnd:ju_bnd+k2d,           &
     &                       kl_bnd:ku_bnd,                             &
     &                       maxblocksf)
      real,save ::  facevarz(nbndvar,                                   &
     &                       il_bnd:iu_bnd,jl_bnd:ju_bnd,               &
     &                       kl_bnd:ku_bnd+k3d,                         &
     &                       maxblocksf)
      integer,save :: interp_mask_facex(nbndvar)
      integer,save :: interp_mask_facey(nbndvar)
      integer,save :: interp_mask_facez(nbndvar)
      integer,save :: interp_mask_facex_res(nbndvar)
      integer,save :: interp_mask_facey_res(nbndvar)
      integer,save :: interp_mask_facez_res(nbndvar)
      integer, save :: gcell_on_fc_pointer(3,nbndvar)
      logical, save :: gcell_on_fc(3,nbndvar)
      logical, save :: int_gcell_on_fc(3,nbndvar)
      logical, save :: checkp_on_fc(3,nbndvar)
#else
      real,allocatable,save ::  facevarx(:,:,:,:,:)
      real,allocatable,save ::  facevary(:,:,:,:,:)
      real,allocatable,save ::  facevarz(:,:,:,:,:)
      integer,allocatable,save :: interp_mask_facex(:)
      integer,allocatable,save :: interp_mask_facey(:)
      integer,allocatable,save :: interp_mask_facez(:)
      integer,allocatable,save :: interp_mask_facex_res(:)
      integer,allocatable,save :: interp_mask_facey_res(:)
      integer,allocatable,save :: interp_mask_facez_res(:)
      integer,allocatable,save :: gcell_on_fc_pointer(:,:)
      logical,allocatable,save :: gcell_on_fc(:,:)
      logical,allocatable,save :: int_gcell_on_fc(:,:)
      logical,allocatable,save :: checkp_on_fc(:,:)
#endif
      integer, save :: ngcell_on_fc(3)
      target :: facevarx,facevary,facevarz


!!****v* physicaldata/facevars
!!
!! NAME
!!
!!   facevarx, facevary, facevarz
!! 
!! SYNOPSIS
!!
!!   public,real ::  facevarx(nbndvar,                                   
!!                            il_bnd:iu_bnd+1,jl_bnd:ju_bnd,             
!!                            kl_bnd:ku_bnd,                             
!!                            maxblocksf)
!!   public,real ::  facevary(nbndvar,                                   
!!                            il_bnd:iu_bnd,jl_bnd:ju_bnd+k2d,           
!!                            kl_bnd:ku_bnd,                             
!!                            maxblocksf)
!!   public,real ::  facevarz(nbndvar,                                   
!!                            il_bnd:iu_bnd,jl_bnd:ju_bnd,               
!!                            kl_bnd:ku_bnd+k3d,                         
!!                            maxblocksf)
!!
!!   public, real, allocatable :: facevarx(:,:,:,:,:)
!!   public, real, allocatable :: facevary(:,:,:,:,:)
!!   public, real, allocatable :: facevarz(:,:,:,:,:)
!!
!! DESCRIPTION
!!
!!   Arrays for holding face centered data.  May or may not have permanent 
!!   guardcell storage allocated. If the LIBRARY flag is defined in 
!!   paramesh_preprocessor.fh these arrays are allocated at runtime.
!!
!!***

!!****v* physicaldata/interp_mask_facex(y,z)
!!
!! NAME
!!
!!   interp_mask_facex, interp_mask_facey, interp_mask_facez
!! 
!! SYNOPSIS
!!
!!   public, integer :: interp_mask_facex(nbndvar)
!!   public, integer :: interp_mask_facey(nbndvar)
!!   public, integer :: interp_mask_facez(nbndvar)
!!
!!   public, integer, allocatable :: interp_mask_facex(:)
!!   public, integer, allocatable :: interp_mask_facey(:)
!!   public, integer, allocatable :: interp_mask_facez(:)
!!
!! DESCRIPTION
!!
!!   Arrays to select different interpolation orders (or other user defined 
!!   interpolation schemes) used during prolongation operations on the face
!!   centered data stored in facevarx(y,z).  The values stored in 
!!   interp_mask_facex(y,z) may be changed at anytime during execution and 
!!   they can take on different values for each different variable stored in 
!!   facevarx(y,z).
!!
!!   The default value(s) for interp_mask_facex(y,z) is 1 for all variables, 
!!   1 through nbndvar.
!!   The default behaviour for using different values of facevarx(y,z) is:
!!     interp_mask_facex(y,z) =  0  -> direct injection
!!     interp_mask_facex(y,z) =  1  -> linear interpolation
!!     interp_mask_facex(y,z) =  2  -> Lagrange polynomial of 2nd order
!!     interp_mask_facex(y,z) =  3  -> Lagrange polynomial of 3rd order
!!     interp_mask_facex(y,z) =  4  -> Lagrange polynomial of 4th order
!!     interp_mask_facex(y,z) >= 20 -> user defined
!!
!!   Usage Example:
!!     If you wish to use the default interpolation during prolongation for 
!!     all variables in facevarx(y,z), but to use a 2nd order polynomial 
!!     interpolation for variable 4 in facevarx then set,
!!
!!     interp_mask_facex(4) = 2
!!
!!***

!!****v* physicaldata/interp_mask_facex(y,z)_res
!!
!! NAME
!!
!!   interp_mask_facex_res, interp_mask_facey_res, interp_mask_facez_res
!! 
!! SYNOPSIS
!!
!!   public, integer :: interp_mask_facex_res(nbndvar)
!!   public, integer :: interp_mask_facey_res(nbndvar)
!!   public, integer :: interp_mask_facey_res(nbndvar)
!!
!!   public, integer, allocatable :: interp_mask_facex_res(:)
!!   public, integer, allocatable :: interp_mask_facey_res(:)
!!   public, integer, allocatable :: interp_mask_facez_res(:)
!!
!! DESCRIPTION
!!
!!   Array to select different interpolation orders (or other user defined 
!!   interpolation schemes) used during restriction operations on the face
!!   centered data stored in facevarx(y,z).  The values stored in 
!!   interp_mask_facex(y,z)_res may be changed at anytime during execution and 
!!   they can take on different values for each different variable stored in 
!!   facevarx(y,z).
!!
!!   The default value(s) for interp_mask_facex(y,z)_res is 1 for all 
!!   variables, 1 through nbndvar.
!!   The default behaviour for using different values of facevarx(y,z) is:
!!     interp_mask_facex(y,z)_res =  0  -> illegal value
!!     interp_mask_facex(y,z)_res =  1  -> Simple average of child data to 
!!                                         parent (linear interpolation)
!!     interp_mask_facex(y,z)_res =  2  -> Lagrange polynomial of 2nd order
!!     interp_mask_facex(y,z)_res =  3  -> Lagrange polynomial of 3rd order
!!     interp_mask_facex(y,z)_res =  4  -> Lagrange polynomial of 4th order
!!     interp_mask_facex(y,z)_res >= 20 -> user defined
!!
!!   Usage Example:
!!     If you wish to use the default interpolation during restriction for 
!!     all variables in facevarx(y,z), but to use a 2nd order polynomial 
!!     interpolation for variable 4 in facevarx then set,
!!
!!     interp_mask_facex(4) = 2
!!
!!***

!!****v* physicaldata/gcell_on_fc
!!
!! NAME
!!
!!   gcell_on_fc
!! 
!! SYNOPSIS
!!
!!   public, logical :: gcell_on_fc(3,nbndvar)
!!   public, logical, allocatable :: gcell_on_fc(:,:)
!!
!! DESCRIPTION
!!
!!   A masking array used to control which variables in the facevarx(y,z) 
!!   array have their guardcells filled during a call to amr_guardcell or 
!!   amr_1blk_guardcell.  If gcell_on_fc(:,ivar) is set to .true. then the 
!!   guardcells for variable 'ivar' in facevarx(y,z) will be filled, if it is 
!!   set to .false. then variable 'ivar' will be skipped and its guardcells 
!!   will not be filled.  
!!
!!   The values in gcell_on_fc can be changed at any point during a run and 
!!   allow a user to avoid unecessary data movement in cases when only a small
!!   set of the variables in facevarx(y,z) need guardcell values at some point
!!   during their algorithm.
!!
!!***

!!****iv physicaldata/int_gcell_on_fc
!!
!! NAME
!!
!!   int_gcell_on_fc
!! 
!! SYNOPSIS
!!
!!   public, logical :: int_gcell_on_fc(3,nbndvar)
!!   public, logical, allocatable :: int_gcell_on_fc(:,:)
!!
!! DESCRIPTION
!!
!!   ??? I don't know ???
!!
!!***

!!****iv physicaldata/ngcell_on_fc
!!
!! NAME
!!
!!   ngcell_on_fc
!! 
!! SYNOPSIS
!!
!!   public, integer :: ngcell_on_fc(3)
!!
!! DESCRIPTION
!!
!!   An integer variable which holds the current number of variables in 
!!   facevarx(y,z) selected for guardcell filling by the mask array 
!!   'gcell_on_fc'.
!!
!!***

!!****iv physicaldata/gcell_on_fc_pointer
!!
!! NAME
!!
!!   gcell_on_fc_pointer
!! 
!! SYNOPSIS
!!
!!   public, integer :: gcell_on_fc_pointer(3,nbndvar)
!!   public, integer, allocatable :: gcell_on_fc_pointer(:,:)
!!
!! DESCRIPTION
!!
!!   I don't know ???
!!
!!
!!***

!!***v* physicaldata/checkp_on_fc
!!
!! NAME
!!
!!   checkp_on_fc
!! 
!! SYNOPSIS
!!
!!   public, logical :: checkp_on_fc(3,nbndvar)
!!   public, logical, allocatable :: checkp_on_fc(:,:)
!!
!! DESCRIPTION
!!
!!   A logical mask array which controls which variables are added to a 
!!   checkpoint by calling the 'amr_checkpoint_wr' routine.  If 
!!   checkp_on_fc(:,ivar) is set to true, then the ivar'th variable in the
!!   facevarx(y,z) array is added to the checkpoint file.  The default 
!!   behaviour has ALL the values in checkp_on_fc set to .true. so that all 
!!   variables in facevarx(y,z) are written to the checkpoint file.
!!
!!***

!----------------------------------------------------------------------------
! EDGE CENTERED DATA
!----------------------------------------------------------------------------

! the solution for cell-edge-centered quantities.
      public :: unk_e_x,unk_e_y,unk_e_z,interp_mask_ec, & 
     &     interp_mask_ec_res
      public :: gcell_on_ec,int_gcell_on_ec
      public :: ngcell_on_ec
      public :: gcell_on_ec_pointer
      public :: checkp_on_ec
#ifndef LIBRARY
      real,save ::  unk_e_x(nbndvare,                                   &
     &                       il_bnd:iu_bnd,jl_bnd:ju_bnd+k2d,           &
     &                       kl_bnd:ku_bnd+k3d,                         &
     &                       maxblocksue)
      real,save ::  unk_e_y(nbndvare,                                   &
     &                       il_bnd:iu_bnd+1,jl_bnd:ju_bnd,             &
     &                       kl_bnd:ku_bnd+k3d,                         &
     &                       maxblocksue)
      real,save ::  unk_e_z(nbndvare,                                   &
     &                       il_bnd:iu_bnd+1,jl_bnd:ju_bnd+k2d,         &
     &                       kl_bnd:ku_bnd,                             &
     &                       maxblocksue)
      integer, save :: interp_mask_ec(nbndvare)
      integer, save :: interp_mask_ec_res(nbndvare)
      integer, save :: gcell_on_ec_pointer(3,nbndvare)
      logical, save :: gcell_on_ec(3,nbndvare)
      logical, save :: int_gcell_on_ec(3,nbndvare)
      logical, save :: checkp_on_ec(3,nbndvare)
#else
      real,allocatable,save ::  unk_e_x(:,:,:,:,:)
      real,allocatable,save ::  unk_e_y(:,:,:,:,:)
      real,allocatable,save ::  unk_e_z(:,:,:,:,:)
      integer,allocatable,save :: interp_mask_ec(:)
      integer,allocatable,save :: interp_mask_ec_res(:)
      integer,allocatable,save :: gcell_on_ec_pointer(:,:)
      logical,allocatable,save :: gcell_on_ec(:,:)
      logical,allocatable,save :: int_gcell_on_ec(:,:)
      logical,allocatable,save :: checkp_on_ec(:,:)
#endif
      integer, save :: ngcell_on_ec(3)
      target :: unk_e_x,unk_e_y,unk_e_z

!!****v* physicaldata/unk_e_x(y,z)
!!
!! NAME
!!
!!   unk_e_x, unk_e_y, unk_e_z
!! 
!! SYNOPSIS
!!
!!   public,real ::  unk_e_x(nbndvare,                                   
!!                           il_bnd:iu_bnd,jl_bnd:ju_bnd+k2d,             
!!                           kl_bnd:ku_bnd+k3d,                             
!!                           maxblocksf)
!!   public,real ::  unk_e_y(nbndvare,                                   
!!                            il_bnd:iu_bnd+1,jl_bnd:ju_bnd,           
!!                            kl_bnd:ku_bnd+k3d,                             
!!                            maxblocksf)
!!   public,real ::  unk_e_z(nbndvare,                                   
!!                           il_bnd:iu_bnd+1,jl_bnd:ju_bnd+k2d,               
!!                           kl_bnd:ku_bnd,                         
!!                           maxblocksf)
!!
!!   public, real, allocatable :: unk_e_x(:,:,:,:,:)
!!   public, real, allocatable :: unk_e_y(:,:,:,:,:)
!!   public, real, allocatable :: unk_e_z(:,:,:,:,:)
!!
!! DESCRIPTION
!!
!!   Arrays for holding edge centered data.  May or may not have permanent 
!!   guardcell storage allocated. If the LIBRARY flag is defined in 
!!   paramesh_preprocessor.fh these arrays are allocated at runtime.
!!
!!***

!!****v* physicaldata/interp_mask_ec(y,z)
!!
!! NAME
!!
!!   interp_mask_ec
!! 
!! SYNOPSIS
!!
!!   public, integer :: interp_mask_ec(nbndvare)
!!
!!   public, integer, allocatable :: interp_mask_ec(:)
!!
!! DESCRIPTION
!!
!!   Arrays to select different interpolation orders (or other user defined 
!!   interpolation schemes) used during prolongation operations on the edge
!!   centered data stored in unk_e_x, y and z.  The values stored in 
!!   interp_mask_ec may be changed at anytime during execution and they 
!!   can take on different values for each different variable stored in 
!!   unk_e_x, y and z.
!!
!!   The default value(s) for interp_mask_ec is 1 for all variables, 
!!   1 through nbndvare.
!!   The default behaviour for using different values of unk_e_x(y,z) is:
!!     interp_mask_ec =  0  -> direct injection
!!     interp_mask_ec =  1  -> linear interpolation
!!     interp_mask_ec =  2  -> Lagrange polynomial of 2nd order
!!     interp_mask_ec =  3  -> Lagrange polynomial of 3rd order
!!     interp_mask_ec =  4  -> Lagrange polynomial of 4th order
!!     interp_mask_ec >= 20 -> user defined
!!
!!   Usage Example:
!!     If you wish to use the default interpolation during prolongation for 
!!     all variables in unk_e_x(y,z), but to use a 2nd order polynomial 
!!     interpolation for variable 4 in unk_e_x, y, and z then set,
!!
!!     interp_mask_ec(4) = 2
!!
!! NOTE:  Using interp_mask_ec controls all 3 arrays unk_e_x, unk_e_y, and
!!        unk_e_z and you do not have the freedom to control the interpolation
!!        order for these arrays individually as you do for the face centered
!!        variables.
!!
!!***

!!****v* physicaldata/interp_mask_ec_res
!!
!! NAME
!!
!!   interp_mask_ec_res
!! 
!! SYNOPSIS
!!
!!   public, integer :: interp_mask_ec_res(nbndvare)
!!
!!   public, integer, allocatable :: interp_mask_ec_res(:)
!!
!! DESCRIPTION
!!
!!   Array to select different interpolation orders (or other user defined 
!!   interpolation schemes) used during restriction operations on the edge
!!   centered data stored in unk_e_x, y, and z.  The values stored in 
!!   interp_mask_ec_res may be changed at anytime during execution and 
!!   they can take on different values for each different variable stored in 
!!   unk_e_x, y, and z.
!!
!!   The default value(s) for interp_mask_ec_res is 1 for all variables, 
!!     1 through nbndvare.
!!   The default behaviour for using different values of unk_e_x(y,z) is:
!!     interp_mask_ec_res =  0  -> illegal value
!!     interp_mask_ec_res =  1  -> Simple average of child data to parent 
!!                                  (linear interpolation)
!!     interp_mask_ec_res =  2  -> Lagrange polynomial of 2nd order
!!     interp_mask_ec_res =  3  -> Lagrange polynomial of 3rd order
!!     interp_mask_ec_res =  4  -> Lagrange polynomial of 4th order
!!     interp_mask_ec_res >= 20 -> user defined
!!
!!   Usage Example:
!!     If you wish to use the default interpolation during restriction for 
!!     all variables in unk_e_x(y,z), but to use a 2nd order polynomial 
!!     interpolation for variable 4 in unk_e_x, y, and z then set,
!!
!!     interp_mask_ec_res(4) = 2
!!
!! NOTE:  Using interp_mask_ec_res controls all 3 arrays unk_e_x, unk_e_y, and
!!        unk_e_z and you do not have the freedom to control the interpolation
!!        order for these arrays individually as you do for the face centered
!!        variables.
!!
!!***

!!****v* physicaldata/gcell_on_ec
!!
!! NAME
!!
!!   gcell_on_ec
!! 
!! SYNOPSIS
!!
!!   public, logical :: gcell_on_ec(3,nbndvare)
!!   public, logical, allocatable :: gcell_on_ec(:,:)
!!
!! DESCRIPTION
!!
!!   A masking array used to control which variables in the unk_e_x(y,z) array
!!   have their guardcells filled during a call to amr_guardcell or 
!!   amr_1blk_guardcell.  If gcell_on_fc(:,ivar) is set to .true. then the 
!!   guardcells for variable 'ivar' in unk_e_x(y,z) will be filled, if it 
!!   is set to .false. then variable 'ivar' will be skipped and its guardcells
!!   will not be filled.  
!!
!!   The values in gcell_on_ec can be changed at any point during a run 
!!   and allow a user to avoid unecessary data movement in cases when only 
!!   a small set of the variables in unk_e_x(y,z) need guardcell values at 
!!   some point during their algorithm.
!!
!!***

!!****iv physicaldata/int_gcell_on_ec
!!
!! NAME
!!
!!   int_gcell_on_ec
!! 
!! SYNOPSIS
!!
!!   public, logical :: int_gcell_on_ec(3,nbndvare)
!!   public, logical, allocatable :: int_gcell_on_ec(:,:)
!!
!! DESCRIPTION
!!
!!   ??? I don't know ???
!!
!!***

!!****iv physicaldata/ngcell_on_ec
!!
!! NAME
!!
!!   ngcell_on_ec
!! 
!! SYNOPSIS
!!
!!   public, integer :: ngcell_on_ec(3)
!!
!! DESCRIPTION
!!
!!   An integer variable which holds the current number of variables in 
!!   unk_e_x(y,z) selected for guardcell filling by the mask array 
!!   'gcell_on_ec'.
!!
!!***

!!****iv physicaldata/gcell_on_ec_pointer
!!
!! NAME
!!
!!   gcell_on_ec_pointer
!! 
!! SYNOPSIS
!!
!!   public, integer :: gcell_on_ec_pointer(3,nbndvare)
!!   public, integer, allocatable :: gcell_on_ec_pointer(:,:)
!!
!! DESCRIPTION
!!
!!   I don't know ???
!!
!!
!!***

!!***v* physicaldata/checkp_on_ec
!!
!! NAME
!!
!!   checkp_on_ec
!! 
!! SYNOPSIS
!!
!!   public, logical :: checkp_on_ec(3,nbndvare)
!!   public, logical, allocatable :: checkp_on_ec(:,:)
!!
!! DESCRIPTION
!!
!!   A logical mask array which controls which variables are added to a 
!!   checkpoint by calling the 'amr_checkpoint_wr' routine.  If 
!!   checkp_on_ec(:,ivar) is set to to true, then the ivar'th variable in 
!!   the unk_e_x(y,z) array is added to the checkpoint file.  The default 
!!   behaviour has ALL the values in checkp_on_ec set to .true. so that all 
!!   variables in unk_e_x(y,z) are written to the checkpoint file.
!!
!!***

!-----------------------------------------------------------------------------
! CELL CORNER DATA
!-----------------------------------------------------------------------------

! the solution for cell-corner based quantities.
      public :: unk_n, interp_mask_nc, interp_mask_nc_res
      public :: gcell_on_nc,int_gcell_on_nc
      public :: ngcell_on_nc
      public :: gcell_on_nc_pointer
      public :: checkp_on_nc
#ifndef LIBRARY
      real,save ::  unk_n(nbndvarc,                                     &
     &                    il_bnd:iu_bnd+1,jl_bnd:ju_bnd+k2d,            &
     &                    kl_bnd:ku_bnd+k3d,                            &
     &                    maxblocksn)
      integer,save :: interp_mask_nc(nbndvarc)
      integer,save :: interp_mask_nc_res(nbndvarc)
      integer, save :: gcell_on_nc_pointer(nbndvarc)
      logical, save :: gcell_on_nc(nbndvarc)
      logical, save :: int_gcell_on_nc(nbndvarc)
      logical, save :: checkp_on_nc(nbndvarc)
#else
      real,allocatable,save ::  unk_n(:,:,:,:,:)
      integer,allocatable,save :: interp_mask_nc(:)
      integer,allocatable,save :: interp_mask_nc_res(:)
      integer,allocatable,save :: gcell_on_nc_pointer(:)
      logical,allocatable,save :: gcell_on_nc(:)
      logical,allocatable,save :: int_gcell_on_nc(:)
      logical,allocatable,save :: checkp_on_nc(:)
#endif

      integer, save :: ngcell_on_nc
      target :: unk_n

!!****v* physicaldata/unk_n
!!
!! NAME
!!
!!   unk_n
!! 
!! SYNOPSIS
!!
!!   public,real ::  unk_n(nbndvarc,                                   
!!                         il_bnd:iu_bnd+1,jl_bnd:ju_bnd+k2d,             
!!                         kl_bnd:ku_bnd+k3d,                             
!!                         maxblocksf)
!!
!!   public, real, allocatable :: unk_n(:,:,:,:,:)
!!
!! DESCRIPTION
!!
!!   Arrays for holding cell corner data.  May or may not have permanent 
!!   guardcell storage allocated. If the LIBRARY flag is defined in 
!!   paramesh_preprocessor.fh these arrays are allocated at runtime.
!!
!!***

!!****v* physicaldata/interp_mask_nc
!!
!! NAME
!!
!!   interp_mask_nc
!! 
!! SYNOPSIS
!!
!!   public, integer :: interp_mask_nc(nbndvarc)
!!
!!   public, integer, allocatable :: interp_mask_nc(:)
!!
!! DESCRIPTION
!!
!!   Arrays to select different interpolation orders (or other user defined 
!!   interpolation schemes) used during prolongation operations on the cell
!!   corner data stored in unk_n.  The values stored in interp_mask_nc may be 
!!   changed at anytime during execution and they can take on different values
!!   for each different variable stored in unk_n.
!!
!!   The default value(s) for interp_mask_nc is 1 for all variables, 
!!   1 through nbndvarc.
!!   The default behaviour for using different values of unk_n is:
!!     interp_mask_nc =  0  -> direct injection
!!     interp_mask_nc =  1  -> linear interpolation
!!     interp_mask_nc =  2  -> Lagrange polynomial of 2nd order
!!     interp_mask_nc =  3  -> Lagrange polynomial of 3rd order
!!     interp_mask_nc =  4  -> Lagrange polynomial of 4th order
!!     interp_mask_nc >= 20 -> user defined
!!
!!   Usage Example:
!!     If you wish to use the default interpolation during prolongation for 
!!     all variables in unk_n, but to use a 2nd order polynomial 
!!     interpolation for variable 4 in unk_n then set,
!!
!!     interp_mask_nc(4) = 2
!!
!!***

!!****v* physicaldata/interp_mask_nc_res
!!
!! NAME
!!
!!   interp_mask_nc_res
!! 
!! SYNOPSIS
!!
!!   public, integer :: interp_mask_nc_res(nbndvarc)
!!
!!   public, integer, allocatable :: interp_mask_nc_res(:)
!!
!! DESCRIPTION
!!
!!   Array to select different interpolation orders (or other user defined 
!!   interpolation schemes) used during restriction operations on the cell
!!   corner data stored in unk_n.  The values stored in interp_mask_nc_res may
!!   be changed at anytime during execution and they can take on different 
!!   values for each different variable stored in unk_n.
!!
!!   The default value(s) for interp_mask_nc_res is 1 for all variables, 
!!     1 through nbndvarc.
!!   The default behaviour for using different values of unk_n is:
!!     interp_mask_nc_res =  0  -> illegal value
!!     interp_mask_nc_res =  1  -> Simple average of child data to parent 
!!                                 (linear interpolation)
!!     interp_mask_nc_res =  2  -> Lagrange polynomial of 2nd order
!!     interp_mask_nc_res =  3  -> Lagrange polynomial of 3rd order
!!     interp_mask_nc_res =  4  -> Lagrange polynomial of 4th order
!!     interp_mask_nc_res >= 20 -> user defined
!!
!!   Usage Example:
!!     If you wish to use the default interpolation during restriction for 
!!     all variables in unk_n, but to use a 2nd order polynomial 
!!     interpolation for variable 4 in unk_n then set,
!!
!!     interp_mask_nc_res(4) = 2
!!
!!
!!***

!!****v* physicaldata/gcell_on_nc
!!
!! NAME
!!
!!   gcell_on_nc
!! 
!! SYNOPSIS
!!
!!   public, logical :: gcell_on_nc(nbndvarc)
!!   public, logical, allocatable :: gcell_on_nc(:)
!!
!! DESCRIPTION
!!
!!   A masking array used to control which variables in the unk_e_n array have 
!!   their guardcells filled during a call to amr_guardcell or 
!!   amr_1blk_guardcell.  If gcell_on_nc(ivar) is set to .true. then the 
!!   guardcells for variable 'ivar' in unk_n will be filled, if it is set to 
!!   .false. then variable 'ivar' will be skipped and its guardcells will not 
!!   be filled.  
!!
!!   The values in gcell_on_nc can be changed at any point during a run and 
!!   allow a user to avoid unecessary data movement in cases when only a small
!!   set of the variables in unk_n need guardcell values at some point during 
!!   their algorithm.
!!
!!***

!!****iv physicaldata/int_gcell_on_nc
!!
!! NAME
!!
!!   int_gcell_on_nc
!! 
!! SYNOPSIS
!!
!!   public, logical :: int_gcell_on_nc(nbndvarc)
!!   public, logical, allocatable :: int_gcell_on_nc(:)
!!
!! DESCRIPTION
!!
!!   ??? I don't know ???
!!
!!***

!!****iv physicaldata/ngcell_on_nc
!!
!! NAME
!!
!!   ngcell_on_nc
!! 
!! SYNOPSIS
!!
!!   public, integer :: ngcell_on_nc(3)
!!
!! DESCRIPTION
!!
!!   An integer variable which holds the current number of variables in unk_n
!!   selected for guardcell filling by the mask array 'gcell_on_nc'.
!!
!!***

!!****iv physicaldata/gcell_on_nc_pointer
!!
!! NAME
!!
!!   gcell_on_nc_pointer
!! 
!! SYNOPSIS
!!
!!   public, integer :: gcell_on_nc_pointer(3,nbndvarc)
!!   public, integer, allocatable :: gcell_on_nc_pointer(:)
!!
!! DESCRIPTION
!!
!!   I don't know ???
!!
!!
!!***

!!***v* physicaldata/checkp_on_nc
!!
!! NAME
!!
!!   checkp_on_nc
!! 
!! SYNOPSIS
!!
!!   public, logical :: checkp_on_nc(nbndvarc)
!!   public, logical, allocatable :: checkp_on_nc(:)
!!
!! DESCRIPTION
!!
!!   A logical mask array which controls which variables are added to a 
!!   checkpoint by calling the 'amr_checkpoint_wr' routine.  If 
!!   checkp_on_nc(ivar) is set to true, then the ivar'th variable in the unk_n
!!   array is added to the checkpoint file.  The default behaviour has ALL the
!!   values in checkp_on_nc set to .true. so that all variables in unk_e_n are
!!   written to the checkpoint file.
!!
!!***

!-----------------
! Timestep control
!-----------------

! arrays used for timestep control
      public :: time_loc,dtlevel,phase_dt,loc_cycle,ncyc_local
      public :: ldtcomplete
      real, save    :: dtlevel(maxlevels)
      integer, save :: phase_dt(maxlevels),loc_cycle(maxlevels)
      integer, save :: ncyc_local(maxlevels)
#ifndef LIBRARY
      logical, save :: ldtcomplete(maxblocks_alloc)
      real, save    :: time_loc(maxblocks_alloc)
#else
      logical, allocatable,save :: ldtcomplete(:)
      real, allocatable,save    :: time_loc(:)
#endif


      public :: t_unk,tfacevarx,tfacevary,tfacevarz
      public :: t_unk_e_x,t_unk_e_y,t_unk_e_z,t_unk_n
!#ifndef LIBRARY
!#if defined(VAR_DT) || defined(PRED_CORR)
!      real, save    :: t_unk(nvar,il_bnd:iu_bnd,jl_bnd:ju_bnd,          &
!     &                            kl_bnd:ku_bnd,                        &
!     &                            maxblocks)
!      real, save  :: tfacevarx(nbndvar,il_bnd:iu_bnd+1,jl_bnd:ju_bnd,   &
!     &                          kl_bnd:ku_bnd,maxblocksf) 
!      real, save  :: tfacevary(nbndvar,il_bnd:iu_bnd,jl_bnd:ju_bnd+k2d, &
!     &                          kl_bnd:ku_bnd,maxblocksf)
!      real, save  :: tfacevarz(nbndvar,il_bnd:iu_bnd,jl_bnd:ju_bnd,     &
!     &                          kl_bnd:ku_bnd+k3d,maxblocksf)
!
!      real, save  :: t_unk_e_x(nbndvare,il_bnd:iu_bnd,jl_bnd:ju_bnd+k2d,&
!     &                          kl_bnd:ku_bnd+k3d,maxblocksue) 
!      real, save  :: t_unk_e_y(nbndvare,il_bnd:iu_bnd+1,jl_bnd:ju_bnd,  &
!     &                          kl_bnd:ku_bnd+k3d,maxblocksue) 
!      real, save  :: t_unk_e_z(nbndvare,il_bnd:iu_bnd+1,                &
!     &                          jl_bnd:ju_bnd+k2d,                      &
!     &                          kl_bnd:ku_bnd,maxblocksue) 
!      real, save  :: t_unk_n(nbndvarc,il_bnd:iu_bnd+1,                  &
!     &                          jl_bnd:ju_bnd+k2d,                      &
!     &                          kl_bnd:ku_bnd+k3d,maxblocksn) 
!#endif
!#else
      real, allocatable,save  :: t_unk(:,:,:,:,:)
      real, allocatable,save  :: tfacevarx(:,:,:,:,:)
      real, allocatable,save  :: tfacevary(:,:,:,:,:)
      real, allocatable,save  :: tfacevarz(:,:,:,:,:)
      real, allocatable,save  :: t_unk_e_x(:,:,:,:,:)
      real, allocatable,save  :: t_unk_e_y(:,:,:,:,:)
      real, allocatable,save  :: t_unk_e_z(:,:,:,:,:)
      real, allocatable,save  :: t_unk_n(:,:,:,:,:)
!#endif /* LIBRARY */


!-----------------------------------------------------------------
! include header file defining 1blk data structure

!------------------------------------------------------------------------------
! data_1blk
!------------------------------------------------------------------------------
!
! This file declares the storage space used to handle the `current
! working block' when the user decides not to reserve permanent
! storage space for guardcells for all blocks, but instead to 
! fill guardcells as needed. This strategy requires 2 working blocks,
! one for the leaf node and one for its parent.

      public :: unk1,facevarx1,facevary1,facevarz1
      public :: unk_e_x1,unk_e_y1,unk_e_z1
      public :: unk_n1

#ifndef LIBRARY
! the solution for cell-centered quantities.
      real, save :: unk1(nvar,il_bnd1:iu_bnd1,jl_bnd1:ju_bnd1,          &
     &                        kl_bnd1:ku_bnd1,                          &
     &                        npblks)
! the solution for cell-face-centered quantities.
      real, save :: facevarx1(nbndvar,il_bnd1:iu_bnd1+1,                &
     &                                jl_bnd1:ju_bnd1,                  &
     &                                kl_bnd1:ku_bnd1,npblks)
      real, save :: facevary1(nbndvar,il_bnd1:iu_bnd1,                  &
     &                                jl_bnd1:ju_bnd1+k2d,              &
     &                                kl_bnd1:ku_bnd1,npblks)
      real, save :: facevarz1(nbndvar,il_bnd1:iu_bnd1,                  &
     &                                jl_bnd1:ju_bnd1,                  &
     &                                kl_bnd1:ku_bnd1+k3d,npblks)  
! the solution for cell-edge-centered quantities.
      real,save ::  unk_e_x1(nbndvare,                                  &
     &                       il_bnd1:iu_bnd1,jl_bnd1:ju_bnd1+k2d,       &
     &                       kl_bnd1:ku_bnd1+k3d,                       &
     &                       npblks)
      real,save ::  unk_e_y1(nbndvare,                                  &
     &                       il_bnd1:iu_bnd1+1,jl_bnd1:ju_bnd1,         &
     &                       kl_bnd1:ku_bnd1+k3d,                        &
     &                       npblks)
      real,save ::  unk_e_z1(nbndvare,                                  &
     &                       il_bnd1:iu_bnd1+1,jl_bnd1:ju_bnd1+k2d,     &
     &                       kl_bnd1:ku_bnd1,                           &
     &                       npblks)
! the solution for cell-corner based quantities.
      real,save ::  unk_n1(nbndvarc,                                    &
     &                       il_bnd1:iu_bnd1+1,jl_bnd1:ju_bnd1+k2d,     &
     &                       kl_bnd1:ku_bnd1+k3d,                       &
     &                       npblks)
#else
! the solution for cell-centered quantities.
      real, save, allocatable :: unk1(:,:,:,:,:)
! the solution for cell-face-centered quantities.
      real, save, allocatable :: facevarx1(:,:,:,:,:)
      real, save, allocatable :: facevary1(:,:,:,:,:)
      real, save, allocatable :: facevarz1(:,:,:,:,:)
! the solution for cell-edge-centered quantities.
      real, save, allocatable ::  unk_e_x1(:,:,:,:,:)
      real, save, allocatable ::  unk_e_y1(:,:,:,:,:)
      real, save, allocatable ::  unk_e_z1(:,:,:,:,:)
! the solution for cell-corner based quantities.
      real, save, allocatable ::  unk_n1(:,:,:,:,:)
#endif

      target :: unk1
      target :: facevarx1, facevary1, facevarz1
      target :: unk_e_x1, unk_e_y1, unk_e_z1
      target :: unk_n1

!------------------------------------------------------------------------------
! workspace arrays for a finer block data layer, to be used in 
! amr_1blk_guardcell_f_to_c
!------------------------------------------------------------------------------

#ifdef LIBRARY
! for cell-centered quantities.
      real, save, allocatable,public ::  unk1_fl(:,:,:,:)
! for cell-face quantities.
      real, save, allocatable,public ::  facevarx1_fl(:,:,:,:)
      real, save, allocatable,public ::  facevary1_fl(:,:,:,:)
      real, save, allocatable,public ::  facevarz1_fl(:,:,:,:)
! for cell-edge quantities.
      real, save, allocatable,public ::  unk_e_x1_fl(:,:,:,:)
      real, save, allocatable,public ::  unk_e_y1_fl(:,:,:,:)
      real, save, allocatable,public ::  unk_e_z1_fl(:,:,:,:)
! for cell-corner quantities.
      real, save, allocatable,public ::  unk_n1_fl(:,:,:,:)
#else
      real, save, public :: unk1_fl(nvar,                               &
     &                              il_bnd1:iu_bnd1+nxb+2*nguard,       &
     &                              jl_bnd1:ju_bnd1+(nyb+2*nguard)*k2d, &
     &                              kl_bnd1:ku_bnd1+(nzb+2*nguard)*k3d )
      real, save, public :: facevarx1_fl(nbndvar,                       &
     &                            il_bnd1:iu_bnd1+nxb+2*nguard+1,       &
     &                            jl_bnd1:ju_bnd1+(nyb+2*nguard)*k2d,   &
     &                            kl_bnd1:ku_bnd1+(nzb+2*nguard)*k3d )
      real, save, public :: facevary1_fl(nbndvar,                       &
     &                            il_bnd1:iu_bnd1+nxb+2*nguard,         &
     &                            jl_bnd1:ju_bnd1+(nyb+2*nguard+1)*k2d, &
     &                            kl_bnd1:ku_bnd1+(nzb+2*nguard)*k3d )
      real, save, public :: facevarz1_fl(nbndvar,                       &
     &                            il_bnd1:iu_bnd1+nxb+2*nguard,         &
     &                            jl_bnd1:ju_bnd1+(nyb+2*nguard)*k2d,   &
     &                            kl_bnd1:ku_bnd1+(nzb+2*nguard+1)*k3d )
      real, save, public :: unk_e_x1_fl(nbndvare,                       &
     &                            il_bnd1:iu_bnd1+nxb+2*nguard,         &
     &                            jl_bnd1:ju_bnd1+(nyb+2*nguard+1)*k2d, &
     &                            kl_bnd1:ku_bnd1+(nzb+2*nguard+1)*k3d )
      real, save, public :: unk_e_y1_fl(nbndvare,                       &
     &                            il_bnd1:iu_bnd1+nxb+2*nguard+1,       &
     &                            jl_bnd1:ju_bnd1+(nyb+2*nguard)*k2d,   &
     &                            kl_bnd1:ku_bnd1+(nzb+2*nguard+1)*k3d )
      real, save, public :: unk_e_z1_fl(nbndvare,                       &
     &                            il_bnd1:iu_bnd1+nxb+2*nguard+1,       &
     &                            jl_bnd1:ju_bnd1+(nyb+2*nguard+1)*k2d, &
     &                            kl_bnd1:ku_bnd1+(nzb+2*nguard)*k3d )
      real, save, public :: unk_n1_fl(nbndvarc, & 
     &                            il_bnd1:iu_bnd1+nxb+2*nguard+1,       &
     &                            jl_bnd1:ju_bnd1+(nyb+2*nguard+1)*k2d, &
     &                            kl_bnd1:ku_bnd1+(nzb+2*nguard+1)*k3d )
#endif

!------------------------------------------------------------------------------
! temporary copy of solution to be used when storing solution prior
! to use of amr_1blk_guardcell
!------------------------------------------------------------------------------

      public :: gt_unk,gt_facevarx,gt_facevary,gt_facevarz
      public :: gt_unk_e_x,gt_unk_e_y,gt_unk_e_z,gt_unk_n
      real, allocatable, save :: gt_unk(:,:,:,:,:)
      real, allocatable, save :: gt_facevarx(:,:,:,:,:)
      real, allocatable, save :: gt_facevary(:,:,:,:,:)
      real, allocatable, save :: gt_facevarz(:,:,:,:,:)
      real, allocatable, save :: gt_unk_e_x(:,:,:,:,:)
      real, allocatable, save :: gt_unk_e_y(:,:,:,:,:)
      real, allocatable, save :: gt_unk_e_z(:,:,:,:,:)
      real, allocatable, save :: gt_unk_n(:,:,:,:,:)
!
! Variables used to control data caching. This helps to avoid unnecessary
! repetition of some guardcell and surrounding-block mapping operations
      public ::  pcache_blk_u,pcache_pe_u,pcache_blk_w,pcache_pe_w
      public ::  lnew_parent
      integer, save :: pcache_blk_u,pcache_pe_u
      integer, save :: pcache_blk_w,pcache_pe_w
      logical, save :: lnew_parent

!-----------------------------------------------------------------
! include header file defining data structure on cell faces

!------------------------------------------------------------------------
! block_boundary_data


!
! This file defines a data structure to be used for quantities
! which may need to be defined at grid block interfaces, eg fluxes,
! pressures.
!

! The convention for relating varaibles associated with cell faces to the
! variables defined at cell centers is as follows:

! If iface_off=0 :
!         the array facevarx(:,i,j,k,:) for example defines data
!         on the x(i-1/2) face of the (i,j,k)-th mesh cell. 
! If iface_off=-1 :
!         the array facevarx(:,i,j,k,:) for example defines data
!         on the x(i+1/2) face of the (i,j,k)-th mesh cell. 


! storage used for fluxes at block boundaries. This is used when conservation
! constraints need to be imposed.

      public :: nfluxvar,nfluxes,maxblocksfl
      public :: flux_x,flux_y,flux_z
      public :: tflux_x,tflux_y,tflux_z
#ifndef LIBRARY
      integer, parameter :: nfluxvar =  N_FLUX_VAR
      integer, parameter :: nfluxes=max(1,nfluxvar)
      integer, parameter :: maxblocksfl=                                &
     &                       1+(maxblocks-1)*min(1,nfluxvar)

      real, save ::  flux_x(nfluxes,1:2,                                &
     &                      jl_bndi:ju_bndi,kl_bndi:ku_bndi,maxblocksfl)
      real, save ::  flux_y(nfluxes,il_bndi:iu_bndi,                      &
     &                              1:2,kl_bndi:ku_bndi,maxblocksfl)
      real, save ::  flux_z(nfluxes,il_bndi:iu_bndi,                      &
     &                              jl_bndi:ju_bndi,1:2,maxblocksfl)
      real, save :: tflux_x(nfluxes,1:2,                                &
     &                              jl_bndi:ju_bndi,                      &
     &                              kl_bndi:ku_bndi,maxblocksfl)
      real, save :: tflux_y(nfluxes,il_bndi:iu_bndi,                      &
     &                              1:2,kl_bndi:ku_bndi,maxblocksfl)
      real, save :: tflux_z(nfluxes,il_bndi:iu_bndi,                      &
     &                              jl_bndi:ju_bndi,1:2,maxblocksfl)

#else
      integer,save :: nfluxvar
      integer,save :: nfluxes
      integer,save :: maxblocksfl
      real, allocatable, save ::  flux_x(:,:,:,:,:)
      real, allocatable, save ::  flux_y(:,:,:,:,:)
      real, allocatable, save ::  flux_z(:,:,:,:,:)
      real, allocatable, save ::  tflux_x(:,:,:,:,:)
      real, allocatable, save ::  tflux_y(:,:,:,:,:)
      real, allocatable, save ::  tflux_z(:,:,:,:,:)
#endif
      target :: flux_x, flux_y, flux_z 


! temporary flux storage needed inside amr_flux_conserve when using
! variable timestep
      public :: ttflux_x,ttflux_y,ttflux_z
      real, allocatable, save :: ttflux_x(:,:,:,:,:)
      real, allocatable, save :: ttflux_y(:,:,:,:,:)
      real, allocatable, save :: ttflux_z(:,:,:,:,:)

! storage used for cell edges at block boundaries. 
! This is used when quantities located at cell edge centers need to
! be used consistently at the boundaries between blocks at different
! refinement levels.

      public :: nedgevar1,nedgevar,nedges,maxblockse
      public :: bedge_facex_y,bedge_facex_z,bedge_facey_x
      public :: bedge_facey_z,bedge_facez_x,bedge_facez_y
      public :: recvarx1e,recvary1e,recvarz1e
      public :: recvarx2e,recvary2e,recvarz2e

#ifndef LIBRARY
      integer, parameter :: nedgevar1=N_EDGE_VAR
      integer, parameter :: nedgevar=max(nedgevar1,nvaredge)
      integer, parameter :: nedges=max(1,nedgevar)

      integer, parameter :: maxblockse=                                 &
     &                           1+(maxblocks-1)*min(1,nedgevar)
      real, save :: bedge_facex_y(nedges,1:2,jl_bnd:ju_bnd+1,           &
     &                                   kl_bnd:ku_bnd+1,maxblockse)
      real, save ::  bedge_facex_z(nedges,1:2,jl_bnd:ju_bnd+1,          &
     &                                    kl_bnd:ku_bnd+1,maxblockse)
      real, save ::  bedge_facey_x(nedges,il_bnd:iu_bnd+1,1:2,          &
     &                                    kl_bnd:ku_bnd+1,maxblockse)
      real, save ::  bedge_facey_z(nedges,il_bnd:iu_bnd+1,1:2,          &
     &                                    kl_bnd:ku_bnd+1,maxblockse)
      real, save ::  bedge_facez_x(nedges,il_bnd:iu_bnd+1,              &
     &                             jl_bnd:ju_bnd+1,1:2,maxblockse)
      real, save ::  bedge_facez_y(nedges,il_bnd:iu_bnd+1,              &
     &                                    jl_bnd:ju_bnd+1,              &
     &                                    1:2,maxblockse)
#else
      integer, save :: nedgevar1
      integer, save :: nedgevar
      integer, save :: nedges
      integer, save :: maxblockse
      real, allocatable, save ::  bedge_facex_y(:,:,:,:,:)
      real, allocatable, save ::  bedge_facex_z(:,:,:,:,:)
      real, allocatable, save ::  bedge_facey_x(:,:,:,:,:)
      real, allocatable, save ::  bedge_facey_z(:,:,:,:,:)
      real, allocatable, save ::  bedge_facez_x(:,:,:,:,:)
      real, allocatable, save ::  bedge_facez_y(:,:,:,:,:)
#endif
#ifndef LIBRARY
      real, save ::  recvarx1e(nedges,1:2,jl_bnd:ju_bnd+1,              &
     &                                kl_bnd:ku_bnd+1)
      real, save ::  recvary1e(nedges,il_bnd:iu_bnd+1,1:2,              &
     &                                kl_bnd:ku_bnd+1)
      real, save ::  recvarz1e(nedges,il_bnd:iu_bnd+1,jl_bnd:ju_bnd+1,  &
     &                                1:2)
      real, save ::  recvarx2e(nedges,1:2,jl_bnd:ju_bnd+1,              &
     &                                kl_bnd:ku_bnd+1)
      real, save ::  recvary2e(nedges,il_bnd:iu_bnd+1,1:2,              &
     &                                kl_bnd:ku_bnd+1)
      real, save ::  recvarz2e(nedges,il_bnd:iu_bnd+1,jl_bnd:ju_bnd+1,  &
     &                                1:2)
#else
      real, allocatable, save ::  recvarx1e(:,:,:,:)
      real, allocatable, save ::  recvary1e(:,:,:,:)
      real, allocatable, save ::  recvarz1e(:,:,:,:)
      real, allocatable, save ::  recvarx2e(:,:,:,:)
      real, allocatable, save ::  recvary2e(:,:,:,:)
      real, allocatable, save ::  recvarz2e(:,:,:,:)
#endif

      public :: tbedge_facex_y,tbedge_facex_z,tbedge_facey_x
      public :: tbedge_facey_z,tbedge_facez_x,tbedge_facez_y      
#ifndef LIBRARY
      real, save :: tbedge_facex_y(nedges,1:2,jl_bnd:ju_bnd+1,          &
     &                                    kl_bnd:ku_bnd+1,maxblockse)
      real, save :: tbedge_facex_z(nedges,1:2,jl_bnd:ju_bnd+1,          &
     &                                    kl_bnd:ku_bnd+1,maxblockse)
      real, save :: tbedge_facey_x(nedges,il_bnd:iu_bnd+1,1:2,          &
     &                                    kl_bnd:ku_bnd+1,maxblockse)
      real, save :: tbedge_facey_z(nedges,il_bnd:iu_bnd+1,1:2,          &
     &                                    kl_bnd:ku_bnd+1,maxblockse)
      real, save :: tbedge_facez_x(nedges,il_bnd:iu_bnd+1,              &
     &                                    jl_bnd:ju_bnd+1,              &
     &                                    1:2,maxblockse)
      real, save :: tbedge_facez_y(nedges,il_bnd:iu_bnd+1,              &
     &                                    jl_bnd:ju_bnd+1,              &
     &                                    1:2,maxblockse)
#else
      real, allocatable, save ::  tbedge_facex_y(:,:,:,:,:)
      real, allocatable, save ::  tbedge_facex_z(:,:,:,:,:)
      real, allocatable, save ::  tbedge_facey_x(:,:,:,:,:)
      real, allocatable, save ::  tbedge_facey_z(:,:,:,:,:)
      real, allocatable, save ::  tbedge_facez_x(:,:,:,:,:)
      real, allocatable, save ::  tbedge_facez_y(:,:,:,:,:)
#endif
      public :: ttbedge_facex_y,ttbedge_facex_z,ttbedge_facey_x
      public :: ttbedge_facey_z,ttbedge_facez_x,ttbedge_facez_y      
      real, allocatable, save ::  ttbedge_facex_y(:,:,:,:,:)
      real, allocatable, save ::  ttbedge_facex_z(:,:,:,:,:)
      real, allocatable, save ::  ttbedge_facey_x(:,:,:,:,:)
      real, allocatable, save ::  ttbedge_facey_z(:,:,:,:,:)
      real, allocatable, save ::  ttbedge_facez_x(:,:,:,:,:)
      real, allocatable, save ::  ttbedge_facez_y(:,:,:,:,:)




! arrays used to store geometry information for the working block
      public :: cell_vol
      public :: cell_area1,cell_area2,cell_area3
      public :: cell_leng1,cell_leng2,cell_leng3
      public :: cell_face_coord1,cell_face_coord2,cell_face_coord3
!#ifndef LIBRARY
!      real :: cell_vol(il_bnd1:iu_bnd1,jl_bnd1:ju_bnd1,kl_bnd1:ku_bnd1)
!      real :: cell_area1(il_bnd1:iu_bnd1+1,jl_bnd1:ju_bnd1,             &
!     &                                        kl_bnd1:ku_bnd1)
!      real :: cell_area2(il_bnd1:iu_bnd1,jl_bnd1:ju_bnd1+k2d,           &
!     &                                        kl_bnd1:ku_bnd1)
!      real :: cell_area3(il_bnd1:iu_bnd1,jl_bnd1:ju_bnd1,               &
!     &                                      kl_bnd1:ku_bnd1+k3d)
!      real :: cell_leng1(il_bnd1:iu_bnd1,jl_bnd1:ju_bnd1+k2d,           &
!     &                                        kl_bnd1:ku_bnd1+k3d)
!      real :: cell_leng2(il_bnd1:iu_bnd1+1,jl_bnd1:ju_bnd1,             &
!     &                                        kl_bnd1:ku_bnd1+k3d)
!      real :: cell_leng3(il_bnd1:iu_bnd1+1,jl_bnd1:ju_bnd1+k2d,         &
!     &                                        kl_bnd1:ku_bnd1)
!      real :: cell_face_coord1(il_bnd1:iu_bnd1+1)
!      real :: cell_face_coord2(jl_bnd1:ju_bnd1+k2d)
!      real :: cell_face_coord3(kl_bnd1:ku_bnd1+k3d)
!#else
      real, allocatable  :: cell_vol(:,:,:)
      real, allocatable  :: cell_area1(:,:,:)
      real, allocatable  :: cell_area2(:,:,:)
      real, allocatable  :: cell_area3(:,:,:)
      real, allocatable  :: cell_leng1(:,:,:)
      real, allocatable  :: cell_leng2(:,:,:)
      real, allocatable  :: cell_leng3(:,:,:)
      real, allocatable  :: cell_face_coord1(:)
      real, allocatable  :: cell_face_coord2(:)
      real, allocatable  :: cell_face_coord3(:)
!#endif

! workspace arrays used for inter-block communications
      public :: nbndmax
      public :: recvarxf,recvaryf,recvarzf
      public :: bndtempx1,bndtempy1,bndtempz1
#ifndef LIBRARY
      integer, parameter :: nbndmax=max(nbndvar,nfluxes)

      real, save :: recvarxf(nfluxes,1:2,jl_bndi:ju_bndi, & 
     &                       kl_bndi:ku_bndi)
      real, save :: recvaryf(nfluxes,il_bndi:iu_bndi,1:2, & 
     &                       kl_bndi:ku_bndi)
      real, save :: recvarzf(nfluxes,il_bndi:iu_bndi,jl_bndi:ju_bndi, & 
     &                       1:2)
      real, save :: bndtempx1(nfluxes,1:2,jl_bndi:ju_bndi, & 
     &                        kl_bndi:ku_bndi)
      real, save :: bndtempy1(nfluxes,il_bndi:iu_bndi,1:2, & 
     &                        kl_bndi:ku_bndi)
      real, save :: bndtempz1(nfluxes,il_bndi:iu_bndi,jl_bndi:ju_bndi, & 
     &                        1:2)

#else
      integer :: nbndmax
      real, allocatable, save :: recvarxf(:,:,:,:)
      real, allocatable, save :: recvaryf(:,:,:,:)
      real, allocatable, save :: recvarzf(:,:,:,:)
      real, allocatable, save :: bndtempx1(:,:,:,:)
      real, allocatable, save :: bndtempy1(:,:,:,:)
      real, allocatable, save :: bndtempz1(:,:,:,:)
#endif



! parameters used in communication calls
      public :: len_block_bndx,len_block_bndy,len_block_bndz
      public :: len_block_ex,len_block_ey,len_block_ez
#ifndef LIBRARY
      integer, parameter :: len_block_bndx= & 
     &  2*(ju_bndi-jl_bndi+1)*(ku_bndi-kl_bndi+1)
      integer, parameter :: len_block_bndy= & 
     &  2*(iu_bndi-il_bndi+1)*(ku_bndi-kl_bndi+1)
      integer, parameter :: len_block_bndz= & 
     &  2*(iu_bndi-il_bndi+1)*(ju_bndi-jl_bndi+1)

      integer, parameter :: len_block_ex=2*(ju_bnd+k2d)*(ku_bnd+k3d)
      integer, parameter :: len_block_ey=2*(iu_bnd+1  )*(ku_bnd+k3d)
      integer, parameter :: len_block_ez=2*(iu_bnd+1  )*(ju_bnd+k2d)
#else
      integer :: len_block_bndx
      integer :: len_block_bndy
      integer :: len_block_bndz
      integer :: len_block_ex
      integer :: len_block_ey
      integer :: len_block_ez
#endif

!-----------------------------------------------------------------
! Array used to store variables which make up any divergence
! free fields
      integer,allocatable,save,public :: i_divf_fc_vars(:,:)

!-----------------------------------------------------------------
! Index arrays used in boundary condition routines.
      integer,public :: bc_index_i(2,3,5)
      integer,public :: bc_index_j(2,3,5)
      integer,public :: bc_index_k(2,3,5)


!-----------------------------------------------------------------
! Logical flags required to signal algorithmic states
      logical, public :: lrestrict_in_progress
      logical, public :: lprolong_in_progress
      logical, public :: lguard_in_progress

!-----------------------------------------------------------------
! Logical flag to indicate if restrictionless guardcell filling
! has been selected
      logical,public :: l_f_to_c

! Index arrays used to record destination data values for fine layer
! neighbor guardcells
      integer,public :: f2c_ind_unk(2,3,27)
      integer,public :: f2c_ind_facex(2,3,27)
      integer,public :: f2c_ind_facey(2,3,27)
      integer,public :: f2c_ind_facez(2,3,27)
      integer,public :: f2c_ind_unkex(2,3,27)
      integer,public :: f2c_ind_unkey(2,3,27)
      integer,public :: f2c_ind_unkez(2,3,27)
      integer,public :: f2c_ind_unkn(2,3,27)

!-----------------------------------------------------------------
! Mpi communication pattern identifier
      integer, public :: mpi_pattern_id

!-----------------------------------------------------------------
! Error trapping and management

! To record whether amr_gsurrounding_blks has been called.
      public :: gsurrblks_set
      integer :: gsurrblks_set

! a counter which can be used to keep track of calls to routines
      public :: instance
      integer :: instance

! Diagonals flag
      public :: diagonals
#ifdef LIBRARY
      logical, save :: diagonals
#else /* LIBRARY */
#ifdef DIAGONALS
      logical, parameter :: diagonals = .true.
#else
      logical, parameter :: diagonals = .false.
#endif
#endif /* LIBRARY */

! error checking flag
      public :: amr_error_checking
#ifdef LIBRARY
      logical, save :: amr_error_checking
#else
#ifdef AMR_ERROR_CHECKING
      logical, parameter :: amr_error_checking = .true.
#else
      logical, parameter :: amr_error_checking = .false.
#endif
#endif /* LIBRARY */

! no_permanent_guardcells flag
      public :: no_permanent_guardcells
#ifdef LIBRARY
      logical, save :: no_permanent_guardcells
#else
#ifdef NO_PERMANENT_GUARDCELLS
      logical, parameter :: no_permanent_guardcells = .true.
#else
      logical, parameter :: no_permanent_guardcells = .false.
#endif
#endif /* LIBRARY */

! advance_all_levels flag
      public :: advance_all_levels
#ifdef LIBRARY
      logical, save :: advance_all_levels
#else
#ifdef ADVANCE_ALL_LEVELS
      logical, parameter :: advance_all_levels = .true.
#else
      logical, parameter :: advance_all_levels = .false.
#endif
#endif /* LIBRARY */

! force_consistency flag
      public :: force_consistency
#ifdef LIBRARY
      logical, save :: force_consistency
#else
#ifdef FORCE_CONSISTENCY_AT_SRL_INTERFACES
      logical, parameter :: force_consistency = .true.
#else
      logical, parameter :: force_consistency = .false.
#endif
#endif /* LIBRARY */

! consv_fluxes and consv_flux_densities flag
      public :: consv_fluxes, consv_flux_densities
#ifdef LIBRARY
      logical, save :: consv_fluxes, consv_flux_densities
#else
#ifdef CONSV_FLUXES
      logical, parameter :: consv_fluxes = .true.
      logical, parameter :: consv_flux_densities = .false.
#else
      logical, parameter :: consv_fluxes = .false.
      logical, parameter :: consv_flux_densities = .true.
#endif
#endif /* LIBRARY */

! edge_value and edge_value_integ flag
      public :: edge_value, edge_value_integ
#ifdef LIBRARY
      logical, save :: edge_value, edge_value_integ
#else
#ifdef EDGE_VALUE
#undef EDGE_VALUE_INTEG
      logical, parameter :: edge_value = .true.
      logical, parameter :: edge_value_integ = .false.
#endif
#ifdef EDGE_VALUE_INTEG
#undef EDGE_VALUE
      logical, parameter :: edge_value = .false.
      logical, parameter :: edge_value_integ = .true.
#endif
#ifndef EDGE_VALUE
#ifndef EDGE_VALUE_INTEG
      logical, parameter :: edge_value = .false.
      logical, parameter :: edge_value_integ = .false.
#endif
#endif
#endif /* LIBRARY */

! var_dt flag
      public :: var_dt
#ifdef LIBRARY
      logical, save :: var_dt
#else
#ifdef VAR_DT
      logical, parameter :: var_dt = .true.
#else
      logical, parameter :: var_dt = .false.
#endif
#endif /* LIBRARY */

! pred_corr flag
      public :: pred_corr
#ifdef LIBRARY
      logical, save :: pred_corr
#else
#ifdef PRED_CORR
      logical, parameter :: pred_corr = .true.
#else
      logical, parameter :: pred_corr = .false.
#endif
#endif /* LIBRARY */

! empty_cells flag
      public :: empty_cells
#ifdef LIBRARY
      logical, save :: empty_cells
#else
#ifdef EMPTY_CELLS
      logical, parameter :: empty_cells = .true.
#else
      logical, parameter :: empty_cells = .false.
#endif
#endif /* LIBRARY */

! conserve flag
      public :: conserve
#ifdef LIBRARY
      logical, save :: conserve
#else
#ifdef CONSERVE
      logical, parameter :: conserve = .true.
#else
      logical, parameter :: conserve = .false.
#endif
#endif /* LIBRARY */

! clean_divb flag
      public :: clean_divb
#ifdef LIBRARY
      logical, save :: clean_divb
#else
#ifdef CLEAN_DIVB
      logical, parameter :: clean_divb = .true.
#else
      logical, parameter :: clean_divb = .false.
#endif
#endif /* LIBRARY */

! divergence_free flag
      public :: divergence_free
#ifdef LIBRARY
      logical, save :: divergence_free
#else
#ifdef DIVERGENCE_FREE
      logical, parameter :: divergence_free = .true.
#else
      logical, parameter :: divergence_free = .false.
#endif
#endif /* LIBRARY */

! curvilinear
      public :: curvilinear, cartesian_pm, cylindrical_pm,  & 
     &          spherical_pm, polar_pm, lsingular_line, & 
     &          curvilinear_conserve
#ifdef LIBRARY
      logical, save :: curvilinear, cartesian_pm, cylindrical_pm,  & 
     &                 spherical_pm, polar_pm, lsingular_line, & 
     &                 curvilinear_conserve
#else
#ifdef CURVILINEAR
      logical, parameter :: curvilinear = .true.
#else
      logical, parameter :: curvilinear = .false.
#endif
#ifdef CURVILINEAR_CONSERVE
      logical, parameter :: curvilinear_conserve = .true.
#else
      logical, parameter :: curvilinear_conserve = .false.
#endif
#ifdef CARTESIAN
      logical, parameter :: cartesian_pm = .true.
#else
      logical, parameter :: cartesian_pm = .false.
#endif
#ifdef CYLINDRICAL
      logical, parameter :: cylindrical_pm = .true.
#else
      logical, parameter :: cylindrical_pm = .false.
#endif
#ifdef SPHERICAL
      logical, parameter :: spherical_pm = .true.
#else
      logical, parameter :: spherical_pm = .false.
#endif
#ifdef POLAR
      logical, parameter :: polar_pm = .true.
#else
      logical, parameter :: polar_pm = .false.
#endif
#ifdef SINGULAR_LINE
      logical, parameter :: lsingular_line = .true.
#else
      logical, parameter :: lsingular_line = .false.
#endif
#endif /* LIBRARY */

!-----------------------------------------------------------------

      end module physicaldata
