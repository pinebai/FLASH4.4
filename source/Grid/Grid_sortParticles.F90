!!****f* source/Grid/Grid_sortParticles
!!
!! NAME
!!  Grid_sortParticles
!!
!! SYNOPSIS
!!
!!  call Grid_sortParticles(real(INOUT)    :: dataBuf(:,:),
!!                          integer(IN)    :: props
!!                          integer(INOUT) :: localCount,
!!                          integer(IN)    :: elementTypes,
!!                          integer(IN)    :: maxCount,
!!                          integer(OUT)   :: elementsPerBlk(:,:),
!!                          integer(IN)    :: attrib1,
!!                 optional,integer(IN)    :: attrib2)
!!                    
!!  
!! DESCRIPTION 
!!  
!!  Sorts the elements by block number. There are two types of blocknumber associated with
!!  elements data structures which have valid values, but are not valid blocknumber in the 
!!  mesh. These are "UNKNOWN" and "NONEXISTENT". The sorter finds out the number of blocks
!!  on the current processor in the mesh, and puts all elements associated with blocknumbers 
!!  in the range of 1 to localNumBlocks in the processor into the appropriate bins. In the bin
!!  for localNumBlocks+1 it puts all elements with block number = UNKNOWN, and in 
!!  localNumBloks+2 it put all elements with block number = NONEXISTENT. For any other block 
!!  number in the BLK_PART_PROP field of any element, the routine aborts. 
!! 
!!
!! ARGUMENTS 
!!
!!  dataBuf : List of elements. It is two dimensional real array, the first dimension
!!              represents each element's properties, and second dimension is index to
!!              elements.
!!
!!  props : number of properties of each element in the dataBuf datastructure
!!
!!  localCount : While coming in it contains the current number of elements mapped to
!!                      this processor. After all the data structure movement, the number
!!                      of local elements might change, and the new value is put back into it.
!!  elementTypes  : Count of different types of elements in the simulation
!!  maxCount : This is parameter determined at runtime, and is the maximum number of local
!!                        elements that a simulation expects to have. All the arrays that hold
!!                        particles in the Particles unit are allocated based on this number.
!!                        
!! elementsPerBlk : integer array containing number of elements on each blk.
!!  attrib1       : the primary property of the elements on which to sort them
!!  attrib2       : if present, then the elements are first sorted on attrib2, and then 
!!                  within each group with a given attrib2 value, on the attrib1.
!!
!!
!! NOTES
!!   currently this routine is called only by io_writeParticles in permanent guardcell mode,
!!   and by the routines that map mesh to particles 
!!   CAUTION : This routine should not be called when " gr_ptSourceBuf" needs to have valid values
!!             since it is used as scratch space by this routine.
!!
!!   The algorithm used for sorting requires (and silently assumes) that
!!     o   if attrib1 or attrib2 is block then for each element is 
!!         either in the valid range 1..MAXBLOCKS,
!!         or has the special value NONEXISTENT; and
!!     o   if the either attribute in particles is  TYPE_PART_PROP 
!!         property for each element, it is in the valid range 1..elementTypes.
!!   Elements with a block value of NONEXISTENT will be effectively dropped (by
!!   potentially reusing their storage locations for other elements, and by not counting them
!!   in the 'elementsPerBlk' output array), but this routine does *not* update the localCount
!!   counter; so the caller better take care of updating its view of the number of elements (probably
!!   from the information in 'elementsPerBlk') after Grid_sortElements returns, if some of the
!!   elements may be NONEXISTENT.
!!
!!***
subroutine Grid_sortParticles(dataBuf,props, localCount,elementTypes,maxCount,&
                              elementsPerBlk,attrib1, attrib2)

  implicit none
#include "Flash.h"
  integer,intent(IN) ::  maxCount,props,elementTypes
  integer, intent(INOUT) :: localCount
  real,dimension(props,maxCount),intent(INOUT) :: dataBuf
  integer,dimension(MAXBLOCKS,elementTypes),intent(OUT) :: elementsPerBlk
  integer, intent(IN) :: attrib1
  integer,optional, intent(IN) :: attrib2

  elementsPerBlk(:,:)=0
  return
end subroutine Grid_sortParticles
