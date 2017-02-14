!!****if* source/flashUtilities/sorting/ut_sortOnProcs
!!
!! NAME
!!  ut_sortOnProcs
!!
!! SYNOPSIS
!!
!!  ut_sortOnProcs(integer (IN) :: count,
!!                 integer (IN) :: props,  
!!                 integer (IN) :: attrib,
!!                 integer (IN) :: numProcs,
!!                 real (INOUT) :: storage, 
!!                 real (INOUT) :: workspace, 
!!                 integer (OUT) :: perProc(:),  
!!                 integer (OUT) :: ifNonZeroProc(:),
!!                 integer (OUT):: nonZeroProcsCount)   
!!                    
!!  
!! DESCRIPTION 
!!  
!!  Given a data structure that is a two dimensional array, where the first dimension
!!  represents properties of individual elements in the data structure, 
!!  and the second dimension represents the count of the elements in the data
!!  structure, this routine finds the processor number associated with each element
!!  and sorts the elements by their processor number.  It assumes that
!!  one of the properties in the first dimension contains the processor number, 
!!  where that index is given by the argument "attrib".  This subroutine also returns
!!  metadata about the sorted data structure.  The array perProc contains the count 
!!  of elements in associated with each processor. The array ifNonZeroProc acts as a
!!  logical array, it has a value 1 if there is at least one element associated with
!!  the corresponding processor and a value 0 otherwise. Finally nonZeroProcsCount
!!  contains the count of procs other than myPE that have a non zero number of
!!  elements associated with them.
!!
!! ARGUMENTS 
!!
!! count - The number of particles in the storage buffer
!! props - count of attributes associated with particle data structure
!! attrib - the index of the property that contains processor number information
!! numProcs - number of processor under consideration
!! storage - buffer containing the particles
!! workspace    - temporary storage for particles while processing
!! perProc      - an array of size equal to number of processors
!!                each entry is the count of particles destined for the 
!!                corresponding processor
!! ifNonZeroProc      - same size array as perProc, here value is 1 if there are 
!!                any particles to be send to the corresponding processor,
!!                otherwise the value is zero
!! nonZeroProcsCount    - The count of the number of processors that will receive
!!                particles sent by myPE
!!
!!
!!***


subroutine ut_sortOnProcs(count, props, attrib, numProcs,&
     storage, workspace, &
     perProc, ifNonZeroProc,nonZeroProcsCount)

#include "constants.h"
#include "Flash.h"

  implicit none

  integer, intent(IN) :: count, props, attrib, numProcs
  real,dimension(props,count),intent(INOUT) ::storage,workspace
  integer,dimension(numProcs),intent(OUT) :: perProc, ifNonZeroProc
  integer, intent(OUT) :: nonZeroProcsCount

  integer :: i,j,k,n

  integer, dimension(numProcs) :: pntr

  perProc=0
  nonZeroProcsCount=0
  ifNonZeroProc=0

  if(count> 0) then
     
     k=1
     n=0
     do i = 1,count
        j=int(storage(attrib,i))
        if(j/=NONEXISTENT) then
           n=n+1
           workspace(1:props,n)=storage(1:props,i)
           j=j+1
           perProc(j)=perProc(j)+1
           if(ifNonZeroProc(j)==0) then
              ifNonZeroProc(j)=1
              nonZeroProcsCount=nonZeroProcsCount+1
           end if
        end if
     end do
     
     pntr(1)=1
     do i = 2,numProcs
        pntr(i)=pntr(i-1)+perProc(i-1)
     end do
     
     do i = 1,n
        j=int(workspace(attrib,i))+1
        k = pntr(j)
        storage(1:props,k)=workspace(1:props,i)
        pntr(j)=pntr(j)+1
     end do
  end if
end subroutine ut_sortOnProcs
