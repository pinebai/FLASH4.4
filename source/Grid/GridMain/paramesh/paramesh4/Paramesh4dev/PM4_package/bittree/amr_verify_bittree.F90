subroutine amr_verify_bittree()
  use bittree, only : amr_identify_block, bittree_block_count
  use paramesh_dimensions, only: ndim
  use Paramesh_comm_data, only: amr_mpi_meshComm
  use tree, only: lnblocks, bsize, coord, lrefine, &
                  grid_xmin, grid_ymin, grid_zmin, &
                  grid_xmax, grid_ymax, grid_zmax
  use iso_c_binding, only: c_bool, c_int
  use Driver_interface, only: Driver_abortFlash
  implicit none
  
#include "Flash_mpi.h"

  integer :: b, locb, lev, proc0, proc1, nprocs, ierr
  integer(c_int) :: nb
  integer :: i
  integer :: lcoord(3)
  logical :: invalid
  
  integer, allocatable :: all_nblock(:), all_recv(:), all_disp(:)
  integer, allocatable :: coord_list(:,:), all_coords(:,:)
  
  call MPI_COMM_RANK(amr_mpi_meshComm, proc0, ierr)
  call MPI_COMM_SIZE(amr_mpi_meshComm, nprocs, ierr)
  
  invalid = .false.
  do b=1, lnblocks
    if(ndim >= 1) lcoord(1) = int((coord(1,b)-grid_xmin)/bsize(1,b))
    if(ndim >= 2) lcoord(2) = int((coord(2,b)-grid_ymin)/bsize(2,b))
    if(ndim >= 3) lcoord(3) = int((coord(3,b)-grid_zmin)/bsize(3,b))
    lev = lrefine(b)
    call amr_identify_block(nprocs,lev,lcoord,proc1,locb)
    if(lrefine(b) /= lev .or. proc0 /= proc1 .or. b /= locb) then
      !if(proc0 .eq. 0) then
      !  if(.not. invalid) print *,'BITTREE IS WRONG!!! proc=', proc0
      !  print *, ' actual/bittree: proc:',proc0,proc1,'locb:',b,locb,'lev:',lrefine(b),lev
      !end if
      invalid = .true.
    end if
  end do
  
  call MPI_Allreduce(MPI_IN_PLACE, invalid, 1, FLASH_LOGICAL, MPI_LOR, amr_mpi_meshComm, ierr)
  
  if(invalid) then
    if(proc0 == 0) allocate(all_nblock(nprocs))
    call MPI_Gather(lnblocks, 1, FLASH_INTEGER, all_nblock, 1, FLASH_INTEGER, 0, amr_mpi_meshComm, ierr)
    if(proc0 == 0) print *, 'gather 1'
    
    allocate(coord_list(1+ndim,lnblocks))
    do b=1, lnblocks
      coord_list(1,b) = lrefine(b)
      if(ndim >= 1) coord_list(1+1,b) = int((coord(1,b)-grid_xmin)/bsize(1,b))
      if(ndim >= 2) coord_list(1+2,b) = int((coord(2,b)-grid_xmin)/bsize(2,b))
      if(ndim >= 3) coord_list(1+3,b) = int((coord(3,b)-grid_xmin)/bsize(3,b))
    end do
    
    if(proc0 == 0) then
      allocate(all_recv(nprocs))
      allocate(all_disp(nprocs))
      allocate(all_coords(1+ndim,sum(all_nblock)))
      all_recv = all_nblock*(1+ndim)
      all_disp(1) = 0
      do i=2, nprocs
        all_disp(i) = all_disp(i-1) + all_recv(i-1)
      end do
    end if
    
    call MPI_Gatherv(coord_list, (1+ndim)*lnblocks, FLASH_INTEGER, all_coords, &
      all_recv, all_disp, FLASH_INTEGER, 0, amr_mpi_meshComm, ierr)
    if(proc0 == 0) print *, 'gather 2'
    
    if(proc0 == 0) then
      all_disp = all_disp/(1+ndim) ! change all_disp to starting 0-based block
      
      open(1349,file='bittree.misery.log',status='NEW')
      
      call bittree_block_count(nb)
      write(1349,*) 'procs=', nprocs, ' actual blocks=', ubound(all_coords,2), ' bittree blocks=', nb
      write(1349,*) ''
      
      do b=1, ubound(all_coords,2)
        if(proc0 < nprocs-1) then
          if(b-1 >= all_disp(proc0+2)) then
            proc0 = proc0 + 1
          end if
        end if
        lev = all_coords(1,b)
        call amr_identify_block(nprocs, lev, all_coords(2:,b), proc1, locb)
        if(lev /= all_coords(1,b) .or. proc0 /= proc1 .or. b-all_disp(proc0+1) /= locb) then
          write(1349,*) 'lev=', all_coords(1,b), ' ijk=', all_coords(2:,b), &
            ' proc=',proc0, ' lblk=', b-all_disp(proc0+1), &
            ' BITTREE lev=',lev,' proc=',proc1,' lblk=', locb
        else
          write(1349,*) 'lev=', all_coords(1,b), ' ijk=', all_coords(2:,b), &
            ' proc=',proc0, ' lblk=', b-all_disp(proc0+1)
        end if
      end do
      close(1349)
      proc0 = 0
    end if
    
    deallocate(coord_list)
    if(proc0 == 0) then
      deallocate(all_nblock)
      deallocate(all_recv)
      deallocate(all_disp)
      deallocate(all_coords)
      call Driver_abortFlash('Bittree suicide, see bittree.misery.log for discrepancies.')
    end if
  end if
end subroutine
