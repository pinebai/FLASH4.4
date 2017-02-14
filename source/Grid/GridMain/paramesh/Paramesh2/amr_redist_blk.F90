subroutine amr_redist_blk(loc,nprocs,mype,dummy)

  ! $RCSfile: amr_redist_blk.F90,v $
  ! $Revision: 1.2 $
  ! $Date: 2004/09/07 00:58:09 $

  use physicaldata, ONLY: nvar, iu_bnd, ju_bnd, ku_bnd, maxblocks, unk
  use Driver_interface, ONLY : Driver_abortFlash
  use tree,         ONLY: maxblocks_tr, newchild, lnblocks,new_lnblocks

  implicit none
  include 'mpif.h'


!arguments
  integer, INTENT(IN)  :: loc(2,maxblocks_tr) 
  integer, INTENT(IN)  :: nprocs, mype        
  integer, INTENT(IN)  :: dummy               

!locals
  integer              :: blocksize             
  integer              :: loc_perm_arr(max(lnblocks,new_lnblocks))

  integer, ALLOCATABLE :: glob_perm_arr(:)            
  integer, ALLOCATABLE :: glob_to_loc_index(:,:)
  integer, ALLOCATABLE :: loc_to_glob_index(:,:)

  integer              :: lnblocks_list(nprocs) 
  integer              :: mystrt(nprocs)

  integer              :: lnblocks_tot, lnblocks_max
  integer              :: i, j, n, count, stat
  integer              :: npos_loc


 

! compute doubles in block size
  blocksize = iu_bnd*ju_bnd*ku_bnd*nvar


! get number of entries in perm matrix for this proc
! careful, arrays may grow or shrink locally on a proc
! after data movement, so we must create a _complete_
! permutation matrix, even for dummy data
  npos_loc = max(new_lnblocks, lnblocks)

! maxblocks check added here. In theory, we should just have to check new_lnblocks,
! but I'll check lnblocks also just in case something went wrong at the previous
! refinement. It can't hurt.
 if (npos_loc .gt. maxblocks) then
     call Driver_abortFlash("fatal error in amr_redist_blk: maxblocks exceeded")
  end if


! send complete list to all procs
  call MPI_Allgather(npos_loc, 1, MPI_INTEGER, lnblocks_list, &
       1, MPI_INTEGER, MPI_COMM_WORLD, stat)


! get total number of blocks across all procs
  lnblocks_tot = sum(lnblocks_list)
  lnblocks_max = maxval(lnblocks_list)

  
! now that we know total blocks, allocate memory
! for global permutation array and global <-> local
! index conversion arrays
  allocate(glob_perm_arr(lnblocks_tot))
  allocate(glob_to_loc_index(2,lnblocks_tot))
  allocate(loc_to_glob_index(lnblocks_max, nprocs))

! list of global to local index converstions
! global index i resides at local index glob_to_loc_index(1,i) 
! on proc glob_to_loc_index(2,i)

  count = 1
  do j = 1, nprocs
     n = lnblocks_list(j)
     do i = 1,n
        glob_to_loc_index(1, count) = i
        glob_to_loc_index(2, count) = j - 1
        loc_to_glob_index(i,j) = count
        count = count + 1
     end do
  end do



! create local version of permutation array
  do i = 1,npos_loc
     if (loc(1,i) .eq. -1) then
        loc_perm_arr(i) = -1
     else
        loc_perm_arr(i) = loc_to_glob_index(loc(1,i),loc(2,i)+1)
     end if
  end do

! get the list of new locs for each proc on all procs
! be careful: lnblocks different on each proc
! store into new arrays pos_glob/proc_glob of size lnblocks_tot
  do n = 1,nprocs
     mystrt(n) = loc_to_glob_index(1,n) - 1
  end do

! get global permutation array by combining locals 
  call MPI_AllgatherV(loc_perm_arr, npos_loc, MPI_INTEGER, glob_perm_arr, &
       lnblocks_list, mystrt, MPI_INTEGER, MPI_COMM_WORLD, stat)

! fill in the blanks 
  call complete_perm_arr(glob_perm_arr,lnblocks_tot)

! do global reorder
  call mpi_reorder(unk, glob_perm_arr, glob_to_loc_index, npos_loc, &
			lnblocks_tot, blocksize, mype)


! deallocate memory
  deallocate(glob_perm_arr)
  deallocate(glob_to_loc_index)
  deallocate(loc_to_glob_index)


end subroutine amr_redist_blk


!------------------------------
!! reorder n blocks of sz doubles
!------------------------------
subroutine mpi_reorder(data, p, glob_to_loc, nloc, n, isz, mype)
  implicit none

  interface
     subroutine pswap(loc_data, iglob, jglob, iloc, jloc, iproc, jproc, isz, mype) 
       real,      dimension(:) :: loc_data
       integer                 :: iglob, jglob, iloc, jloc, iproc, jproc, isz, mype
     end subroutine pswap
  end interface

  integer,   intent(IN)                      :: nloc, n, isz, mype
  real,    intent(INOUT), dimension(nloc*isz):: data
  integer, intent(INOUT), dimension(n)       :: p
  integer, intent(IN), dimension(2,n)        :: glob_to_loc


  integer :: i, j, k

!! inverse of p
  integer,   dimension(n) :: q


!! calculate p inverse
  do i = 1, n
     q(p(i)) = i
  end do


!! main reorder loop
  do i = 1,n
     j = p(i)
     k = q(i)
     call pswap(data, i, k, glob_to_loc(1,i), glob_to_loc(1,k),  &
          &     glob_to_loc(2,i), glob_to_loc(2,k), isz, mype)
     p(k) = p(i)
     p(i) = i
     q(j) = q(i)
     q(i) = i
  end do
    
end subroutine mpi_reorder

!-------------------------------------
!! swap two locs in same address space
!-------------------------------------
subroutine swap(x, i, j, isz, tmp) 
  real, dimension(:)      :: x
  integer                 :: i, j, isz
  real, dimension(isz)    :: tmp
 
!! locals
  integer                 :: istrt, iend, jstrt, jend


  istrt = (i-1)*isz   + 1
  iend  = istrt + isz - 1
  jstrt = (j-1)*isz   + 1
  jend  = jstrt + isz - 1

  tmp = x(istrt:iend)
  x(istrt:iend) = x(jstrt:jend)
  x(jstrt:jend) = tmp

end subroutine swap

!-------------------------------------
!! swap two blobs in possibly different 
!! address spaces 
!-------------------------------------
subroutine pswap(loc_data, iglob, jglob, iloc, jloc, iproc, jproc, isz, mype) 
  implicit none
  include 'mpif.h'

  interface
     subroutine swap(x, i, j, isz, tmp) 
       real, dimension(:)      :: x
       integer                 :: i, j, isz
       real, dimension(isz)    :: tmp
     end subroutine swap
  end interface

  real,      dimension(:) :: loc_data
  integer                 :: iglob, jglob, iloc, jloc, iproc, jproc, isz, mype
  
!! one block worth of scratch space
  real, dimension(isz)    :: tmp

  integer                 :: istrt, iend, jstrt, jend
  
  integer                 :: ierr, status(MPI_STATUS_SIZE)
  

!! if global to and from indices are identical do nothing 
  if (iglob .eq. jglob) return	
  
!! if global to and from indices reside on same PE, do local swap 
  if (iproc .eq. jproc .and. jproc .eq. mype) then
     call swap(loc_data, iloc, jloc, isz, tmp)
     return
  end if
  
!! if I am proc participting in one end of swap, send data residing in
!! corresponding local index -- iloc -- to destination proc -- jproc.
!! At the same time, receive a message from jproc and store in _tmp.

  istrt = (iloc-1)*isz   + 1
  iend  = istrt + isz - 1

  jstrt = (jloc-1)*isz   + 1
  jend  = jstrt + isz - 1


  if (iproc .eq. mype) then
     call MPI_SENDRECV(loc_data(istrt), isz, MPI_DOUBLE_PRECISION, jproc, iloc,    &
		      tmp,             isz, MPI_DOUBLE_PRECISION, jproc, jloc,    &
		      MPI_COMM_WORLD, status, ierr)
     loc_data(istrt:iend) = tmp
  end if


!! if I am proc participting in one end of swap, send data residing in
!! corresponding local index -- iloc -- to destination proc -- jproc.
!! At the same time, receive a message from jproc and store in _tmp.

  if (jproc .eq. mype) then
    call MPI_SENDRECV(loc_data(jstrt), isz, MPI_DOUBLE_PRECISION, iproc, jloc, &
		      tmp,             isz, MPI_DOUBLE_PRECISION, iproc, iloc, &
                      MPI_COMM_WORLD, status, ierr)
    loc_data(jstrt:jend) = tmp
 end if

  return
end subroutine pswap

!-------------------------------------
!! naive algorithm to fill in -1's in glob_perm array
!! with non-represented indices so we have true complete
!! permuation array, which mpi_sort requires
!! This is O[n^2] but the numbers are tiny so no real
!! need to optimize
!-------------------------------------

subroutine complete_perm_arr(glob_perm_arr, n)
  implicit none
  integer, intent(IN) :: n
  integer, dimension(n), intent(INOUT) :: glob_perm_arr
  
  integer :: i,j,next
  logical :: hit
  integer, dimension(n) :: missing

  next = 0


  do i = 1,n
     hit = .false.
     do j = 1,n
        if (glob_perm_arr(j) .eq. i) then
           hit = .true.
           exit
        end if
     end do
     if (.not. hit) then
        next = next + 1
        missing(next) = i
     end if
  end do

  next = 0
  do i = 1,n
     if (glob_perm_arr(i) .eq. -1) then
        next = next + 1
        glob_perm_arr(i) = missing(next)
     end if
  end do

end subroutine complete_perm_arr
