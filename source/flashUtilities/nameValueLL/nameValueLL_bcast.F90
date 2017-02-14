!!****if* source/flashUtilities/nameValueLL/nameValueLL_bcast
!!
!! NAME
!!  nameValueLL_bcast
!!
!! SYNOPSIS
!!  
!!  nameValueLL_bcast (context_type(INOUT) :: context,
!!                     integer(IN)         :: myPE)
!!  
!! DESCRIPTION 
!!  
!!  Broadcasts any context from a namedValueLinkedList from the
!!  MASTER_PE to the other processors.  Broadcast the parameter names
!!  and values to all other processors For each parameter, first
!!  broadcast its type separately, then use the appropriate
!!  single-parameter broadcast routine to send the parameter
!!  information.
!!
!! ARGUMENTS
!!
!!   context -- structure holding all the data to be broadcast
!!   myPE    -- current processor number
!!
!!***

subroutine nameValueLL_bcast(context, myPE)

  use nameValueLL_data !, ONLY: context_type, nameValueLL_add, &
!     &  name_invalid, name_real, name_int, name_str, name_log, &
!     &  real_list_type, int_list_type, str_list_type, log_list_type, TYPE_VAR
  use Driver_interface, ONLY : Driver_abortFlash

  implicit none

#include "constants.h"
#include "Flash_mpi.h"  

  type (context_type), intent(inout)          :: context    
  integer, intent(in)                         :: myPE

  type (real_list_type), pointer:: node_real
  type (int_list_type), pointer :: node_int
  type (str_list_type), pointer :: node_str
  type (log_list_type), pointer :: node_log
  integer                           :: istat, i, strtype
  integer                           :: n_real, n_int, n_str, n_log 
  real, allocatable                 :: real_vals(:)
  integer, allocatable              :: int_vals(:)
  logical, allocatable              :: log_vals(:)
  character(len=MAX_STRING_LENGTH), allocatable    :: str_vals(:), & 
       &                                       real_names(:), & 
       &                                       int_names(:), & 
       &                                       str_names(:), & 
       &                                       log_names(:)


  ! Create an MPI derived datatype to send character strings.

  call MPI_Type_Contiguous (MAX_STRING_LENGTH, MPI_CHARACTER, strtype, istat)
  call MPI_Type_Commit (strtype, istat)



  ! Only the MASTER_PE processor has the parameters; it sends all of
  ! the parameters of each type in turn to the rest of the
  ! processors.


  ! Send the real-valued parameters.

  n_real = context%n_real
  call MPI_Bcast (n_real, 1, MPI_INTEGER, MASTER_PE, & 
       &                  MPI_COMM_WORLD, istat)

  if (n_real > 0) then

     allocate (real_names(n_real), real_vals(n_real), stat=istat)
     if (istat /= 0) then
        write (*,*) 'nameValueLL_bcast :  allocate failed in real bcast'
        call Driver_abortFlash("Error: nameValueLL_bcast :  allocate failed in real bcast")            
     endif

     if (myPE == MASTER_PE) then
        i = 1
        node_real => context%real_list
        do while (associated(node_real))
           real_names(i) = node_real%name
           real_vals(i)  = node_real%value
           i = i + 1
           node_real => node_real%next
        enddo
     endif

     call MPI_Bcast (real_names, n_real, & 
          &                    strtype, MASTER_PE, MPI_COMM_WORLD, istat)
     call MPI_Bcast (real_vals, n_real, & 
          &                    FLASH_REAL, MASTER_PE, & 
          &                    MPI_COMM_WORLD, istat)

     if (myPE /= MASTER_PE) then
        do i = 1, n_real
           call nameValueLL_add(context, real_names(i), real_vals(i), TYPE_VAR)
        enddo
     endif

     deallocate (real_names, real_vals)

  endif


  ! Send the integer-valued parameters.

  n_int = context%n_int
  call MPI_Bcast (n_int, 1, MPI_INTEGER, MASTER_PE, & 
       &                  MPI_COMM_WORLD, istat)

  if (n_int > 0) then

     allocate (int_names(n_int), int_vals(n_int), stat=istat)
     if (istat /= 0) then
        write (*,*) 'nameValueLL_bcast:  allocate failed'
        call Driver_abortFlash("Error: nameValueLL_bcast :  allocate failed")
     endif

     if (myPE == MASTER_PE) then
        i = 1
        node_int => context%int_list
        do while (associated(node_int))
           int_names(i) = node_int%name
           int_vals(i)  = node_int%value
           i = i + 1
           node_int => node_int%next
        enddo
     endif

     call MPI_Bcast (int_names, n_int, & 
          &                    strtype, MASTER_PE, MPI_COMM_WORLD, istat)
     call MPI_Bcast (int_vals, n_int, & 
          &                    MPI_INTEGER, MASTER_PE, & 
          &                    MPI_COMM_WORLD, istat)

     if (myPE /= MASTER_PE) then
        do i = 1, n_int
           call nameValueLL_add (context, int_names(i), int_vals(i), TYPE_VAR)
        enddo
     endif

     deallocate (int_names, int_vals)

  endif

  !  Send the string-valued parameters.

  n_str = context%n_str
  call MPI_Bcast (n_str, 1, MPI_INTEGER, MASTER_PE, & 
       &                  MPI_COMM_WORLD, istat)

  if (n_str > 0) then

     allocate (str_names(n_str), str_vals(n_str), stat=istat)
     if (istat /= 0) then
        write (*,*) 'nameValueLL_bcast :  allocate failed'
        call Driver_abortFlash("Error: nameValueLL_bcast :  allocate failed");
     endif

     if (myPE == MASTER_PE) then
        i = 1
        node_str => context%str_list
        do while (associated(node_str))
           str_names(i) = node_str%name
           str_vals(i)  = node_str%value
           i = i + 1
           node_str => node_str%next
        enddo
     endif

     call MPI_Bcast (str_names, n_str, & 
          &                    strtype, MASTER_PE, MPI_COMM_WORLD, istat)
     call MPI_Bcast (str_vals, n_str, & 
          &                    strtype, MASTER_PE, MPI_COMM_WORLD, istat)

     if (myPE /= MASTER_PE) then
        do i = 1, n_str
           call nameValueLL_add (context, str_names(i), str_vals(i), TYPE_VAR)
        enddo
     endif

     deallocate (str_names, str_vals)

  endif

  !  Send the logical-valued parameters.

  n_log = context%n_log
  call MPI_Bcast (n_log, 1, MPI_INTEGER, MASTER_PE, & 
       &                  MPI_COMM_WORLD, istat)

  if (n_log > 0) then

     allocate (log_names(n_log), log_vals(n_log), stat=istat)
     if (istat /= 0) then
        write (*,*) 'nameValueLL_bcast :  allocate failed'
        call Driver_abortFlash("nameValueLL_bcast :  allocate failed")
     endif

     if (myPE == MASTER_PE) then
        i = 1
        node_log => context%log_list
        do while (associated(node_log))
           log_names(i) = node_log%name
           log_vals(i)  = node_log%value
           i = i + 1
           node_log => node_log%next
        enddo
     endif

     call MPI_Bcast (log_names, n_log, & 
          &                    strtype, MASTER_PE, MPI_COMM_WORLD, istat)
     call MPI_Bcast (log_vals, n_log, & 
          &                    MPI_LOGICAL, MASTER_PE, & 
          &                    MPI_COMM_WORLD, istat)

     if (myPE /= MASTER_PE) then
        do i = 1, n_log
           call nameValueLL_add(context, log_names(i), log_vals(i), TYPE_VAR)
        enddo
     endif

     deallocate (log_names, log_vals)

  endif

  !  Free up the MPI derived datatype used to send strings.

  call MPI_Type_Free (strtype, istat)

  return
end subroutine nameValueLL_bcast
