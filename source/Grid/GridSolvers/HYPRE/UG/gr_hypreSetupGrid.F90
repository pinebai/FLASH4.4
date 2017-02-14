!!****if* source/Grid/GridSolvers/HYPRE/UG/gr_hypreSetupGrid
!!
!!  NAME
!!
!!   gr_hypreSetupGrid
!!
!!  SYNOPSIS
!!
!!  call gr_hypreSetupGrid (integer(IN)  :: blockCount,
!!                          integer(IN), dimension(blockCount) :: blockList,
!!                 OPTIONAL,integer(IN)  :: nvars)
!!
!!  DESCRIPTION 
!! This routine sets up the HYPRE Grid, Called only once in UG. 
!! 
!!
!! ARGUMENTS
!!   blockCount     : The number of blocks in the list.   
!!   blockList      : The list of blocks on which the solution must be updated.
!!   nvars          : Number of variables, also number of equations, for a
!!                    system. Default is 1.
!!
!! SIDE EFFECTS
!!
!!  
!! NOTES
!!
!!   Uses HYPRE library.  
!!
!! SEE ALSO
!!
!!  gr_hypreGridStatus
!!
!!***

subroutine gr_hypreSetupGrid (blockCount, blockList, nvars)

  use Grid_interface, ONLY : Grid_getBlkIndexLimits, Grid_getBlkCornerID  

  use gr_hypreData,     ONLY : gr_hypreGridIsSetUp, gr_hypreLower, gr_hypreUpper, gr_hypreVecX, &
                               gr_hypreVecB, gr_hypreMatA, gr_hypreGraph, gr_hypreStencils, gr_hypreNParts
  use gr_hypreData,     ONLY : gr_hypreNVars, gr_hypreGrid, gr_hypreSolverType,gr_hypreRefineMIN, gr_hypreRefineMAX,&
                               gr_hypreAnisoDiffusion

  use Grid_data,        ONLY : gr_meshComm, gr_meshMe
  use Driver_interface, ONLY : Driver_abortFlash
  
  implicit none
  
#include "constants.h"   
#include "HYPREf.h"  
#include "Flash.h" 
  
  integer,                      intent(IN) :: blockCount
  integer,dimension(blockCount),intent(IN) :: blockList
  integer,OPTIONAL,             intent(IN) :: nvars
  
!!$  integer :: offsets(2*NDIM+1,NDIM)  
  integer, allocatable :: offsets(:,:)
  integer :: datasize(MDIM)
  integer :: lb, blockID 
  integer :: ilower(MDIM), iupper(MDIM)
  integer :: stride (MDIM), lb_stride(blockCount, MDIM)
  integer :: part
  integer :: object_type
  integer :: ent, var, i
  integer :: nentries
  integer :: ierr
  integer, dimension(2,MDIM) :: blkLimitsGC, blkLimits
  
  !!     This comes from 'sstruct_mv/HYPRE_sstruct_mv.h'
  integer ::   HYPRE_SSTRUCT_VARIABLE_CELL = 0 !! CELL CENTERED   
  
  !! Do nothing, if grid already setup.  
  if (gr_hypreGridIsSetUp) return  
  
  gr_hypreGridIsSetUp = .TRUE.
  
  if (gr_hypreSolverType == HYPRE_SPLIT) then     
     object_type = HYPRE_SSTRUCT
  else
     object_type = HYPRE_PARCSR
  end if
  
  !! Lower and upper Corner ID of each local leaf block
  allocate (gr_hypreLower(blockCount, NDIM))
  allocate (gr_hypreUpper(blockCount, NDIM))      
  
  !!-----------------------------------------------------------------------
  !!     1.  Store lower corner ID, upper corner ID of local leaf blocks
  !!     2.  Compute max/min refinement levels on all leaf blocks.
  !!     3.  Store stride of all local leaf blocks.        
  !!-----------------------------------------------------------------------
  do lb=1, blockCount
     blockID = blockList(lb)  
     call Grid_getBlkCornerID(blockId, ilower(1:MDIM), stride)    
     call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)         
     ilower(1:MDIM) =  ceiling(real(ilower(1:MDIM)) / real(stride(1:MDIM)))          
     datasize  (1:MDIM)= blkLimits(HIGH,1:MDIM)-blkLimits(LOW,1:MDIM) 
     iupper(1:MDIM) = ilower(1:MDIM) + datasize (1:MDIM)         
     
     do i=1, NDIM
        gr_hypreLower(lb,i) =  ilower(NDIM-i+1)
        gr_hypreUpper(lb,i) =  iupper(NDIM-i+1)
     end do
     
     lb_stride(lb,1:MDIM) = stride(1:MDIM)       
     
  end do
  
  gr_hypreRefineMIN = 1
  gr_hypreRefineMAX = 1  
  
  !!-----------------------------------------------------------------------
  !!     5.  All blocks (box) at a refinement level forms a part.
  !!         NOTE: Total parts in UG would be 1.
  !!     6.  Only one variable to solve.
  !!-----------------------------------------------------------------------
  gr_hypreNParts = 1
  if (present(nvars)) then
     gr_hypreNVars  = nvars
  else
     gr_hypreNVars  = 1
  end if

  part   = 0   !! part iterator


  !!-----------------------------------------------------------------------
  !!     7.  Create a HYPRE grid object with computed number of 
  !!         parts and given dimension.
  !!-----------------------------------------------------------------------  
  
  call HYPRE_SStructGridCreate(gr_meshComm, NDIM, gr_hypreNParts, gr_hypreGrid, ierr) 
  

  !!-----------------------------------------------------------------------
  !!     8.  Each leaf block in FLASH is a box object (defined by its
  !!         extents) in HYPRE. Add the leaf blocks to their respective parts
  !!         (refinement level in FLASH).
  !!-----------------------------------------------------------------------  
  do lb=1, blockCount 
     call HYPRE_SStructGridSetExtents(gr_hypreGrid, part,gr_hypreLower(lb,1:NDIM),  &
                                      gr_hypreUpper(lb,1:NDIM),ierr)
  end do
  
  !!-----------------------------------------------------------------------
  !!     9.  We are solving for cell centered variables, 
  !!         let HYPRE know that. Also, set number of variables per part. It
  !!         is possible to have more then one variable per part with
  !!         different variable types.
  !!-----------------------------------------------------------------------  
  
  call HYPRE_SStructGridSetVariables(gr_hypreGrid, part, gr_hypreNVars, &
       (/(HYPRE_SSTRUCT_VARIABLE_CELL, var=1,gr_hypreNVars)/) ,ierr)  
  
  !!-----------------------------------------------------------------------
  !!     11. Assemble HYPRE grid (Global call within HYPRE).
  !!-----------------------------------------------------------------------
  call HYPRE_SStructGridAssemble(gr_hypreGrid, ierr)
  
  !!-----------------------------------------------------------------------
  !!     12. Define the discretization stencil, we use standard central
  !!         differencing for diffusion operator, hence we get 3,5,7 point
  !!         stencil for 1D, 2D and 3D respectively (in the isotropic case).
  !!----------------------------------------------------------------------    
  if (.not. gr_hypreAnisoDiffusion) then
     nentries = NDIM*2 + 1
  else
     ! HYPRE for anisotropic diffusion only supports 2D and 3D.
     ! This is meant to be used with multidimensional MHD.
     if (NDIM == 1) then
        call Driver_abortFlash&
          ('The HYPRE GridSolver for anisotropic diffusion is only available for multidimensional setups')
     elseif (NDIM == 2) then
        nentries = 9
     elseif (NDIM == 3) then
        nentries = 19 !(27-8)
     endif
  endif

 
  !!-----------------------------------------------------------------------
  !!     13. Offsets provide relative positions in Matrix, this saves us 
  !!         the effort of computing exact positions on global matrix, 
  !!         transformation of offsets to Matrix is done implicitly by HYPRE. 
  !!         NOTE:  These offsets are different from the order in the HYPRE
  !!                manual: because we use FORTRAN arrays, the indices
  !!                are reversed.
  !!----------------------------------------------------------------------    

  allocate(offsets(nentries,NDIM))
  
#if NDIM == 1
  offsets(1,1) =  0
  offsets(2,1) = -1
  offsets(3,1) =  1
#endif
  
#if NDIM == 2

  if (.not. gr_hypreAnisoDiffusion) then

     !! Isotropic diffusion (isotropic thermal conduction, radiation diffusion, etc.)
     offsets(1,1) =  0
     offsets(1,2) =  0 

     offsets(2,1) =  0
     offsets(2,2) = -1

     offsets(3,1) =  0
     offsets(3,2) =  1

     offsets(4,1) = -1
     offsets(4,2) =  0

     offsets(5,1) =  1
     offsets(5,2) =  0

  else
     !! Anisotropic diffusion; used for
     !! anisotropic thermal conduction along B-field lines in MHD
     offsets(1,1) =  0
     offsets(1,2) =  0

     offsets(2,1) = -1
     offsets(2,2) =  0

     offsets(3,1) =  1
     offsets(3,2) =  0

     offsets(4,1) =  0
     offsets(4,2) = -1

     offsets(5,1) = -1
     offsets(5,2) = -1

     offsets(6,1) =  1
     offsets(6,2) = -1

     offsets(7,1) =  0
     offsets(7,2) =  1

     offsets(8,1) = -1
     offsets(8,2) =  1

     offsets(9,1) =  1
     offsets(9,2) =  1

  endif

#endif
  
#if NDIM == 3
  if (.not. gr_hypreAnisoDiffusion) then

     !! Isotropic diffusion
     offsets(1,1) =  0
     offsets(1,2) =  0 
     offsets(1,3) =  0 

     offsets(2,1) =  0
     offsets(2,2) =  0
     offsets(2,3) = -1

     offsets(3,1) =  0
     offsets(3,2) =  0
     offsets(3,3) =  1
   
     offsets(4,1) =  0
     offsets(4,2) = -1
     offsets(4,3) =  0

     offsets(5,1) =  0
     offsets(5,2) =  1
     offsets(5,3) =  0

     offsets(6,1) = -1
     offsets(6,2) =  0
     offsets(6,3) =  0

     offsets(7,1) =  1
     offsets(7,2) =  0
     offsets(7,3) =  0

  else

     !! Anisotropic diffusion
     offsets(1,1)  =  0
     offsets(1,2)  =  0
     offsets(1,3)  =  0

     offsets(2,1)  = -1
     offsets(2,2)  =  0
     offsets(2,3)  =  0

     offsets(3,1)  =  1
     offsets(3,2)  =  0
     offsets(3,3)  =  0

     offsets(4,1)  =  0
     offsets(4,2)  = -1
     offsets(4,3)  =  0

     offsets(5,1)  = -1
     offsets(5,2)  = -1
     offsets(5,3)  =  0

     offsets(6,1)  =  1
     offsets(6,2)  = -1
     offsets(6,3)  =  0

     offsets(7,1)  =  0
     offsets(7,2)  =  1
     offsets(7,3)  =  0

     offsets(8,1)  = -1
     offsets(8,2)  =  1 
     offsets(8,3)  =  0

     offsets(9,1)  =  1
     offsets(9,2)  =  1
     offsets(9,3)  =  0

     offsets(10,1) =  0
     offsets(10,2) =  0
     offsets(10,3) = -1

     offsets(11,1) = -1
     offsets(11,2) =  0
     offsets(11,3) = -1

     offsets(12,1) =  1
     offsets(12,2) =  0
     offsets(12,3) = -1

     offsets(13,1) =  0
     offsets(13,2) = -1
     offsets(13,3) = -1

     offsets(14,1) =  0
     offsets(14,2) =  1
     offsets(14,3) = -1

     offsets(15,1) =  0
     offsets(15,2) =  0
     offsets(15,3) =  1

     offsets(16,1) = -1
     offsets(16,2) =  0
     offsets(16,3) =  1

     offsets(17,1) =  1
     offsets(17,2) =  0
     offsets(17,3) =  1

     offsets(18,1) =  0
     offsets(18,2) = -1
     offsets(18,3) =  1

     offsets(19,1) =  0
     offsets(19,2) =  1
     offsets(19,3) =  1

  endif
#endif
  
  !!-----------------------------------------------------------------------
  !!     14. Step done for pure convinience, by assigning specific numerical
  !!         values to offsets, we can easily refer to them later (as when
  !!         needed). In the single-variable case, var is just 0 for our problem.
  !!----------------------------------------------------------------------

  allocate(gr_hypreStencils(0:0))
  call HYPRE_SStructStencilCreate(NDIM, nentries, gr_hypreStencils(0), ierr)

  var    = 0   !! var iterator.
  do ent = 1, nentries
     call HYPRE_SStructStencilSetEntry(gr_hypreStencils(0), ent-1, offsets(ent,1:NDIM), var, ierr)
  enddo

  deallocate(offsets)

  !!-----------------------------------------------------------------------
  !!     15. Create Graph object, a graph object is a catch all, it is used
  !!         for building up relations across parts, i.e., at fine-coarse
  !!         boundaries, these cannot be handled using normal stenciled objects.
  !!         NOTE: Not needed in UG.
  !!----------------------------------------------------------------------  
  
  call HYPRE_SStructGraphCreate(gr_meshComm, gr_hypreGrid, gr_hypreGraph, ierr)
  
  
  !!-----------------------------------------------------------------------
  !!     16. HYPRE provides a variety of storage formats. Specific formats
  !!         support specific solvers. HYPRE_PARCSR supports most of the 
  !!         solvers, especially PCG/AMG which we are interested in using.
  !!----------------------------------------------------------------------  
  
  call HYPRE_SStructGraphSetObjectType(gr_hypreGraph,object_type, ierr)  
  
  !!-----------------------------------------------------------------------
  !!     17. Associate a stencil with each part/variable. We use the same
  !!         stencil for all parts.
  !!----------------------------------------------------------------------

  call HYPRE_SStructGraphSetStencil(gr_hypreGraph, part, var, gr_hypreStencils(0), ierr)
  
  !!-----------------------------------------------------------------------
  !!     18. Creating graphs, done if and only if we are using PARAMESH
  !!         and gr_hypreNParts > 1. Which ensure we have a fine-coasrse
  !!         interface.
  !!----------------------------------------------------------------------
 
  !!-----------------------------------------------------------------------
  !!     21. Assemble the HYPRE graph object
  !!-----------------------------------------------------------------------
  call HYPRE_SStructGraphAssemble(gr_hypreGraph, ierr)  
  
  !!-----------------------------------------------------------------------
  !!     22. Create empty matrix & vector objects.
  !!-----------------------------------------------------------------------
  call HYPRE_SStructMatrixCreate(gr_meshComm, gr_hypreGraph, gr_hypreMatA, ierr)      
  call HYPRE_SStructVectorCreate(gr_meshComm, gr_hypreGrid,  gr_hypreVecB, ierr)
  call HYPRE_SStructVectorCreate(gr_meshComm, gr_hypreGrid,  gr_hypreVecX, ierr)  
  
!!$  call HYPRE_SStructMatrixSetSymmetric(gr_hypreMatA, part, var, var, .true., ierr)

  
  !!-----------------------------------------------------------------------
  !!     23. Set storage format type.
  !!         As stated before, this format supports PCG/AMG.
  !!-----------------------------------------------------------------------
  call HYPRE_SStructMatrixSetObjectTyp(gr_hypreMatA, object_type, ierr)
  call HYPRE_SStructVectorSetObjectTyp(gr_hypreVecB, object_type, ierr)
  call HYPRE_SStructVectorSetObjectTyp(gr_hypreVecX, object_type, ierr)
  
  
  !!-----------------------------------------------------------------------
  !!     24. HYPRE allocates memory for objects.
  !!-----------------------------------------------------------------------
  call HYPRE_SStructMatrixInitialize(gr_hypreMatA, ierr) 
  call HYPRE_SStructVectorInitialize(gr_hypreVecB, ierr)
  call HYPRE_SStructVectorInitialize(gr_hypreVecX, ierr)        
  
end subroutine gr_hypreSetupGrid
