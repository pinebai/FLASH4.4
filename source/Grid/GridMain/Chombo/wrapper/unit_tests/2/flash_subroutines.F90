#include "constants.h"
#include "Flash.h"

subroutine driver_init_flash()
  use flash_interfaces, ONLY : Grid_getBlkPtr, Grid_getBlkIndexLimits
  use iso_c_binding, ONLY : c_int, c_double, c_char, c_bool
  use flash_ftypes, ONLY : flash_ug_info_t
  use chombo_f_c_interface, ONLY : ch_define_uniform_grid, ch_fill_guardcells
  use Driver_interface, ONLY : Driver_abortFlash
  implicit none
  include "mpif.h"
  integer (c_int), dimension(MDIM) :: procGrid
  integer (c_int) :: blkID = 1, dim = NDIM, nVar = NUNK_VARS, gds = CENTER, strLen=4
  real (c_double), dimension(MDIM) :: lowDomain, highDomain
  integer :: ierr, myPE, numProcs, v, i, j, k
  integer, dimension(2,MDIM) :: blkLimits,blkLimitsGC
  real, dimension(:,:,:,:), pointer :: dataPtr

  type(flash_ug_info_t) :: flashUGInfo
  integer, parameter :: Total_strings = &
       NUNK_VARS + &
       NFACE_VARS + &
       NFACE_VARS + &
       NFACE_VARS + &
       NSCRATCH_GRID_VARS + &
       NSCRATCH_CENTER_VARS + &
       NSCRATCH_FACEX_VARS + &
       NSCRATCH_FACEY_VARS + &
       NSCRATCH_FACEZ_VARS
  integer, parameter, dimension(MAX_GRID_DATA_STRUCT_TMP) :: Mesh_types = &
       (/ CENTER, &
       FACEX, &
       FACEY, &
       FACEZ, &
       SCRATCH, & 
       SCRATCH_CTR, & 
       SCRATCH_FACEX, &
       SCRATCH_FACEY, &
       SCRATCH_FACEZ /)
  integer, parameter, dimension(MAX_GRID_DATA_STRUCT_TMP) :: Mesh_variables = &
       (/ NUNK_VARS, &
       NFACE_VARS, &
       NFACE_VARS, &
       NFACE_VARS, &
       NSCRATCH_GRID_VARS, &
       NSCRATCH_CENTER_VARS, &
       NSCRATCH_FACEX_VARS, &
       NSCRATCH_FACEY_VARS, &
       NSCRATCH_FACEZ_VARS /)
  integer(kind=c_int), dimension(Total_strings) :: meshStringLen
  character(kind=c_char), dimension &
       (Total_strings * MAX_STRING_LENGTH) :: meshStrings
  logical(kind=c_bool) :: restart
  restart = .false.

  call MPI_Init(ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD, myPE, ierr)
  call MPI_Comm_size(MPI_COMM_WORLD, numProcs, ierr)

  if (numProcs /= NUM_PROCS) then
     call Driver_abortFlash("This is a 4 processor test")
  end if


  !Create a uniform grid distributed across the IAXIS.
  meshStrings(1)='d'; meshStrings(2)='e'; meshStrings(3)='n'; meshStrings(4)='s'
  meshStrings(5)='p'; meshStrings(6)='r'; meshStrings(7)='e'; meshStrings(8)='s'
  meshStringLen = (/strLen,strLen/)
  procGrid = (/numProcs,1,1/)

  flashUGInfo % lowDomain = (/0.0,0.0,0.0/)
  flashUGInfo % highDomain = (/1.0,1.0,1.0/)
  flashUGInfo % procGrid = procGrid
  flashUGInfo % baseDomainSize = (/procGrid(1)*NXB,procGrid(2)*NYB,procGrid(3)*NZB/)
  flashUGInfo % guardCells = (/NGUARD,NGUARD,NGUARD/)
  flashUGInfo % domainBC = (/PERIODIC,PERIODIC,0/)
  flashUGInfo % meshTypes(1:MAX_GRID_DATA_STRUCT_TMP) = &
       Mesh_Types(1:MAX_GRID_DATA_STRUCT_TMP)
  flashUGInfo % meshNumVars(1:MAX_GRID_DATA_STRUCT_TMP) = &
       Mesh_variables(1:MAX_GRID_DATA_STRUCT_TMP)
  flashUGInfo % verbosity = 9

  call ch_define_uniform_grid(flashUGInfo, meshStringLen, meshStrings, restart)


  call Grid_getBlkIndexLimits(blkId, blkLimits, blkLimitsGC)
  call Grid_getBlkPtr(blkID, dataPtr)
  do v = 1, nVar
     do k = blkLimits(LOW,KAXIS), blkLimits(HIGH,KAXIS)
        do j = blkLimits(LOW,JAXIS), blkLimits(HIGH,JAXIS)
           do i = blkLimits(LOW,IAXIS), blkLimits(HIGH,IAXIS)
              dataPtr(i,j,k,v) = real(((v-1)*numProcs)+myPE)
           end do
        end do
     end do
  end do
  !Generally in FLASH we should then release the block pointer

  call ch_fill_guardcells()

end subroutine driver_init_flash


subroutine driver_evolve_flash()
  use flash_interfaces, ONLY : Grid_getBlkPtr, Grid_getBlkIndexLimits
  implicit none
  include "mpif.h"
  real, dimension(:,:,:,:), pointer :: dataPtr
  integer :: blockID, i, j, k, v, myPE, numProcs, ierr
  logical :: success
  integer, dimension(2,MDIM) :: blkLimits,blkLimitsGC
  integer, parameter :: STR_LEN = 200
  character (len=*), parameter :: fileNamePrefix = "Output"
  character (len=STR_LEN) :: fileName, uniqueName, formatStr, formatStrGC
  integer :: iAxisOffset

  !v1 and v2 contain the global expected values (inc. guardcells)
  !for variables (components in Chombo terminology) 1 and 2.
  !We use periodic boundary conditions.  Global arrays are
  !hard-coded for 2D, 4 processors and IAXIS decmposition.
  integer, dimension(NUM_PROCS*(NXB+2*NGUARD),NYB+2*NGUARD) :: v1 = &
       reshape(SOURCE = (/&
       3,0,0,1, 0,1,1,2, 1,2,2,3, 2,3,3,0, &
       3,0,0,1, 0,1,1,2, 1,2,2,3, 2,3,3,0, &
       3,0,0,1, 0,1,1,2, 1,2,2,3, 2,3,3,0, &
       3,0,0,1, 0,1,1,2, 1,2,2,3, 2,3,3,0, &
       3,0,0,1, 0,1,1,2, 1,2,2,3, 2,3,3,0 &
       /), SHAPE = (/NUM_PROCS*(NXB+2*NGUARD),NYB+2*NGUARD/))

  integer, dimension(NUM_PROCS*(NXB+2*NGUARD),NYB+2*NGUARD) :: v2 = &
       reshape(SOURCE = (/&
       7,4,4,5, 4,5,5,6, 5,6,6,7, 6,7,7,4, &
       7,4,4,5, 4,5,5,6, 5,6,6,7, 6,7,7,4, &
       7,4,4,5, 4,5,5,6, 5,6,6,7, 6,7,7,4, &
       7,4,4,5, 4,5,5,6, 5,6,6,7, 6,7,7,4, &
       7,4,4,5, 4,5,5,6, 5,6,6,7, 6,7,7,4 &
       /), SHAPE = (/NUM_PROCS*(NXB+2*NGUARD),NYB+2*NGUARD/))

  !vExpected contain the local processor's portion of v1 and v2.
  integer, dimension(NXB+2*NGUARD,NYB+2*NGUARD) :: vExpected

#ifdef DUMP_LOCAL_GRID
  logical :: dumpLocalGrid = .true.
#else
  logical :: dumpLocalGrid = .false.
#endif


  success = .true.
  blockID = 1;

  call MPI_Comm_rank(MPI_COMM_WORLD, myPE, ierr)
  call MPI_Comm_size(MPI_COMM_WORLD, numProcs, ierr)
  
  call Grid_getBlkIndexLimits(blockId, blkLimits, blkLimitsGC)
  call Grid_getBlkPtr(blockID, dataPtr)


  !Check block limits contain expected values:
  if (blkLimitsGC(HIGH,IAXIS) /= (NXB + (K1D*2*NGUARD))) success = .false.
  if (blkLimitsGC(HIGH,JAXIS) /= (NYB + (K2D*2*NGUARD))) success = .false.
  if (blkLimitsGC(HIGH,KAXIS) /= (NZB + (K3D*2*NGUARD))) success = .false.

  !Check block cells contain expected values:
  iAxisOffset = myPE * (NXB+2*NGUARD)
  do v = lbound(dataPtr,4), ubound(dataPtr,4)
     if (v == 1) then
        vExpected(:,:) = v1(iAxisOffset+1:iAxisOffset+(NXB+2*NGUARD),:)
     else
        vExpected(:,:) = v2(iAxisOffset+1:iAxisOffset+(NXB+2*NGUARD),:)
     end if

     do k = blkLimitsGC(LOW,KAXIS), blkLimitsGC(HIGH,KAXIS)
        do j = blkLimitsGC(LOW,JAXIS), blkLimitsGC(HIGH,JAXIS)
           do i = blkLimitsGC(LOW,IAXIS), blkLimitsGC(HIGH,IAXIS)
              if (int(dataPtr(i,j,k,v)) /= vExpected(i,j)) then
                 success = .false.
              end if
           end do
        end do
     end do
  end do
  !Generally in FLASH we should then release the block pointer

  if (success .eqv. .true.) then
     write(6,'(a,i2,a)') "SUCCESS!  PE ", myPE, &
          " picked up expected C++ assigned values"
  else
     write(6,'(a,i2,a)') "FAILURE!  PE ", myPE, &
          " did not pick up expected C++ assigned values"
  end if



  !Set dumpLocalGird to true if you want to generate a separate
  !file for each processor that shows the local grid points.
  if (dumpLocalGrid .eqv. .true.) then

     !Create format strings to write a row of internal cells and
     !a row of all cells, respectively.
     write (formatStr, fmt='(''('',i3.3,''i3)'')') &
          (blkLimits(HIGH,IAXIS) - blkLimits(LOW,IAXIS) + 1)
     write (formatStrGC, fmt='(''('',i3.3,''i3)'')') &
          (blkLimitsGC(HIGH,IAXIS) - blkLimitsGC(LOW,IAXIS) + 1)


     do v = lbound(dataPtr,4), ubound(dataPtr,4)

        !Create a unique name composed of processor ID and variable.
        !It will look something like: OutputV2_P3.dat.
        write(uniqueName, fmt='(''V'',i1.1,''_P'',i1.1,''.dat'')') v, myPE
        fileName = fileNamePrefix // uniqueName

        open(unit=12, file=fileName)
        write(12,*) "(Local grid) Internal cells:"
        do k = blkLimits(LOW,KAXIS), blkLimits(HIGH,KAXIS)
           do j = blkLimits(LOW,JAXIS), blkLimits(HIGH,JAXIS)
              write(12,formatStr) &
                   int(dataPtr(blkLimits(LOW,IAXIS):blkLimits(HIGH,IAXIS),j,k,v))
           end do
        end do

        write(12,*) "(Local grid) All cells:"
        do k = blkLimitsGC(LOW,KAXIS), blkLimitsGC(HIGH,KAXIS)
           do j = blkLimitsGC(LOW,JAXIS), blkLimitsGC(HIGH,JAXIS)
              write(12,formatStrGC) &
                   int(dataPtr(blkLimitsGC(LOW,IAXIS):blkLimitsGC(HIGH,IAXIS),j,k,v))
           end do
        end do

        if (v == 1) then
           vExpected(:,:) = v1(iAxisOffset+1:iAxisOffset+(NXB+2*NGUARD),:)
        else
           vExpected(:,:) = v2(iAxisOffset+1:iAxisOffset+(NXB+2*NGUARD),:)
        end if

        write(12,*) "(Local grid expected) All cells:"
        do k = blkLimitsGC(LOW,KAXIS), blkLimitsGC(HIGH,KAXIS)
           do j = blkLimitsGC(LOW,JAXIS), blkLimitsGC(HIGH,JAXIS)
              write(12,formatStrGC) &
                   int(vExpected(blkLimitsGC(LOW,IAXIS):blkLimitsGC(HIGH,IAXIS),j))
           end do
        end do
        
        close(12)
     end do
  end if

end subroutine driver_evolve_flash


subroutine driver_finalize_flash()
  use chombo_f_c_interface, ONLY : ch_finalize
  implicit none
  integer :: ierr

  call ch_finalize
  call MPI_Finalize(ierr)
end subroutine driver_finalize_flash
