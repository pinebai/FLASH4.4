!!****if* source/Simulation/SimulationMain/magnetoHD/unitTest/NohCylindricalRagelike/Noh_unitTest
!!
!! NAME
!!
!!  Noh_unitTest
!!
!! SYNOPSIS
!!
!!  use Noh_unitTest, ONLY: Noh_compareDensity
!!
!!  call Noh_unitTest(real(IN)::time, LOGICAL(INOUT)::checkedDensity)
!!
!! DESCRIPTION
!!
!!  This routine compares the output of the MHD Noh
!!  test in cylindrical coordinates to the analytic
!!  solution with a linearly interpolated goodness 
!!  of fit
!!
!! NOTES
!!
!!
!!
!!***

module Noh_unitTest

  contains

    Subroutine Noh_compareDensity(time,checkedDensity)
           use Grid_interface, ONLY : Grid_getBlkIndexLimits, &
                                      Grid_getCellCoords, &
                                      Grid_getBlkPtr, &
                                      Grid_getBlkData, &
                                      Grid_getListOfBlocks, &
                                      Grid_releaseBlkPtr
          use Driver_data, ONLY: dr_globalMe
          use Logfile_interface,   ONLY : Logfile_stamp, Logfile_close
          
#include "constants.h"
#include "Flash.h"

    implicit none
        
    real, intent(IN)                  :: time
    logical, intent(INOUT)            :: checkedDensity
    integer                           :: blockCount
    integer                           :: temp
    real                              :: ChiSquare, interp_slope, interp_b
    real                              :: expect_dens
    integer, dimension(MAXBLOCKS)     :: blockList
    integer, dimension(LOW:HIGH,MDIM) :: blkLimits
    integer, dimension(LOW:HIGH,MDIM) :: blkLimitsGC
    integer                           :: i,n,m,blockID
    integer                           :: ix,iy,iz
    real, pointer, dimension(:,:,:,:) :: U
    real, allocatable                 :: xcent(:)
    real, dimension(0:56)             :: an_dens, an_x
    character(len=20)                 :: fileName
    integer,dimension(4)              :: prNum

! Double check that we are at least past t = .29 and we haven't performed
! this check already
   if ((time .gt. .29) .and. (.not. checkedDensity)) then

    !     Open file of expected values and read to array
    OPEN (1, FILE = "magnoh-analytic.txt", STATUS="OLD") 
    do i=0,56
        read(1,*,err=99) an_x(i),an_dens(i)          
    enddo
    !     Create file names for test result per proc
    temp = dr_globalMe
    do i = 1,4
        prNum(i)= mod(temp,10)
        temp = temp/10
    end do
    filename = "unitTest_"//char(48+prNum(4))//char(48+prNum(3))//&
                                    char(48+prNum(2))//char(48+prNum(1))
    
    open(2,file=fileName)
    write(2,'("P",I0)') dr_globalMe
    
    call Grid_getListOfBlocks(LEAF,blockList,blockCount)


    ChiSquare = 0.0
    !     Begin loop through blocks
    do i=1,blockCount
    
        blockID = blockList(i)
        !     Get block data from current block 
        call Grid_getBlkPtr(blockID,U,CENTER)  
        call Grid_getBlkIndexLimits(blockId,blkLimits,blkLimitsGC)
        allocate(xcent(blkLimits(LOW, IAXIS) : blkLimits(HIGH, IAXIS)))
        call Grid_getCellCoords(IAXIS, blockID, CENTER, .FALSE., &
                    xcent, blkLimits(HIGH, IAXIS)-blkLimits(LOW, IAXIS)+1)
        !     Begin loops through axis to retrieve densities
        do ix = blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS)
          do iy = blkLimits(LOW,JAXIS),blkLimits(LOW,JAXIS)
            do iz = blkLimits(LOW,KAXIS),blkLimits(HIGH,KAXIS)
            
            m=0
            !     Test for cutoff at a distance of .6
            if (xcent(ix) .lt. .6) then
            
                !     Sort the data into the expected array to find closest data
                !     point.  m = index of x value in expected array <= numerical x
                do n = 0,55
                    if ((xcent(ix).ge.an_x(n).and.xcent(ix).lt.an_x(n+1))) m=n
                enddo
                if (xcent(ix).ge.an_x(56)) m=55
                
                ! find interpolated slope of interval around xcent
                interp_slope = (an_dens(m)-an_dens(m+1))/(an_x(m)-an_x(m+1))
                
                ! find y-intercept of linear interpolation of interval
                interp_b = an_dens(m) - an_x(m)*interp_slope 
                
                ! calculate expected dens from linear interpolation
                expect_dens = xcent(ix)*interp_slope + interp_b
                
                ! add error to ChiSquare
                ChiSquare = ChiSquare + (U(DENS_VAR, ix, iy, iz)/100. &
                    - expect_dens)**2 /  expect_dens
                    
                !print *,(U(DENS_VAR,ix,iy,iz)/100-expect_dens)**2/expect_dens, &
                !    ChiSquare   
            endif
            enddo    
          enddo
        enddo
        !     End axis loop
        deallocate(xcent)
        call Grid_releaseBlkPtr(blockID,U,CENTER)
    
        
    enddo
    !     Output to unitTest_#### files whether pass or fail
    if (ChiSquare .le. .16) then
        write(2,'("Noh unitTest PASSED!  all results conformed with expected values.")')
        write(*,'("Noh unitTest PASSED!  all results conformed with expected values.")')
        call Logfile_stamp( "Noh unitTest PASSED!")
        write(2,*) '"ChiSquare" for this process: ', ChiSquare
        if (dr_globalMe==MASTER_PE) &
             print *, '"ChiSquare" for the master process: ', ChiSquare
    else
        write(2,'("Noh unitTest FAILED!")')
        write(*,'("Noh unitTest FAILED!")')
        call Logfile_stamp( "Noh unitTest FAILED!")
        write(2,*) '"ChiSquare" for this process: ', ChiSquare
        print *,   '"ChiSquare" for process',dr_globalMe,': ', ChiSquare
    endif
    ! Claim that this test was already performed
    checkedDensity = .true.
   endif
   goto 98
99 print *, "Reference file magnoh-analytic.txt does not exist"
   close(1)
   close(2)
   call Logfile_close()
98 return
end subroutine Noh_compareDensity
end module
