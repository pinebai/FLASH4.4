!!****if* source/Grid/GridMain/Chombo/AMR/gr_sanitizeDataAfterInterp
!!
!! NAME
!!
!!  gr_sanitizeDataAfterInterp
!!
!!
!! SYNOPSIS
!!
!!  gr_sanitizeDataAfterInterp(integer(in) :: blkList(count),
!!                             integer(in) :: count,
!!                             integer(in) :: layers(MDIM))
!!
!!
!!
!! DESCRIPTION
!!
!!  Given a list of blocks of data, loop over all of the blocks and
!!  check whether solution data in certain variables lie in a reasonable
!!  range of values.
!!
!!  Energies (ENER_VAR and EINT_VAR) are expected to be .ge. gr_smalle,
!!  and the density (DENS_VAR) is expected to be .ge. gr_smallrho,
!!  where gr_smalle and gr_smallrho are lower bounds coming from the
!!  runtime parameters smalle and smlrho, respectively.
!!
!!  For data that do satisfy these expectations, warning messages are
!!  generated, but the offending data is not modified.
!!
!! ARGUMENTS
!! 
!!   blkList - integer list of blocks to be operated on
!!
!!   count - number of blocks in the blkList
!!
!!   layers - number of guardcell layers to be included in the check 
!!
!! NOTES
!!
!!  The checks are based on gr_conserveToPrimitive, which is called to
!!  convert solution data back from conserved form when using the old
!!  (convertToConsvdForMeshCalls) way of ensuring that the mesh
!!  handles data interpolation in conserved form.
!!
!!  This is meant to be called where gr_conserveToPrimitive used to be
!!  called when using the new (convertToConsvdInMeshInterp) way of
!!  ensuring that the mesh handles data interpolation in conserved
!!  form.
!!
!! SEE ALSO
!!
!!  gr_conserveToPrimitive
!!
!!***

!!REORDER(4): solnData

#define DEBUG_CONSCONV

subroutine gr_sanitizeDataAfterInterp(blkList,count, info, layers)

  use Grid_data, ONLY : gr_smallrho, gr_smalle, gr_meshMe
  use Grid_interface, ONLY : Grid_getBlkPtr, Grid_releaseBlkPtr, &
       Grid_getBlkIndexLimits
  use Logfile_interface, ONLY : Logfile_stamp

  implicit none
#include "constants.h"
#undef REAL_FORMAT
#define REAL_FORMAT "(1PG23.16)"
#include "Flash.h"

  integer,intent(IN) :: count
  integer, dimension(count), intent(IN) :: blkList
  character(len=*), intent(IN) :: info
  integer,dimension(MDIM), intent(IN):: layers

  real, dimension(:,:,:,:), pointer :: solnData
  integer,dimension(2,MDIM) :: blkLimits, blkLimitsGC
  integer :: n, blockID
  integer :: myPe, i,j
  integer :: iskip, jskip, kskip
  integer :: il,iu,jl,ju,kl,ku
  character(len=32), dimension(4,2) :: block_buff
  character(len=32)                 :: number_to_str

111 format (a,a,a1,(1x,a18,'=',a),(1x,a2,'=',a5),(1x,a5,'=',a),(1x,a4,'=',a))
112 format (i3,1x,16(1x,1G8.2))


  iskip = NGUARD - layers(IAXIS)
  jskip = (NGUARD - layers(JAXIS)) * K2D
  kskip = (NGUARD - layers(KAXIS)) * K3D


  do n = 1,count

     blockID = blkList(n)
     call Grid_getBlkPtr(blockID, solnData)
     call Grid_getBlkIndexLimits(blockID, blkLimits, blkLimitsGC)

     il = blkLimitsGC(LOW,IAXIS) + iskip
     iu = blkLimitsGC(HIGH,IAXIS) - iskip
     jl = blkLimitsGC(LOW,JAXIS) + jskip
     ju = blkLimitsGC(HIGH,JAXIS) - jskip
     kl = blkLimitsGC(LOW,KAXIS) + kskip
     ku = blkLimitsGC(HIGH,KAXIS) - kskip


#ifdef DENS_VAR
     ! small limits -- in case the interpolants are not monotonic
     if (any(solnData(DENS_VAR,il:iu,jl:ju,kl:ku) .LT. gr_smallrho)) then
        write (block_buff(1,1), '(a18)') 'min. solnData(DENS_VAR)'
        !        write (number_to_str, '('//REAL_FORMAT//',a1)') minval(solnData(DENS_VAR,il:iu,jl:ju,kl:ku,block)), ','
        write (number_to_str, '(G30.22)') minval(solnData(DENS_VAR,il:iu,jl:ju,kl:ku))
        write (block_buff(1,2), '(a)') trim(adjustl(number_to_str))

        write (block_buff(2,1), '(a)') 'PE'
        write (number_to_str, '(i7)') myPe
        write (block_buff(2,2), '(a)') trim(adjustl(number_to_str))

        write (block_buff(3,1), '(a)') 'block'
        write (number_to_str, '(i7)') blockID
        write (block_buff(3,2), '(a)') trim(adjustl(number_to_str))

        write (block_buff(4,1), '(a)') 'type'
        write (number_to_str, '(i7)') LEAF
        write (block_buff(4,2), '(a)') trim(adjustl(number_to_str))

        call Logfile_stamp(block_buff, 4, 2, 'WARNING '//info)
        print 111, 'WARNING ',info,':', ((block_buff(i,j),j=1,2),i=1,4)
#ifdef DEBUG_CONSCONV
        ! For 2D, this prints a slice at the lowest k index that is interior - KW
        do j=blkLimitsGC(HIGH,JAXIS), blkLimitsGC(LOW,JAXIS),-1
           print 112, j, (solnData(DENS_VAR,i,j,blkLimits(LOW,KAXIS)), i=blkLimitsGC(LOW,IAXIS),blkLimitsGC(HIGH,IAXIS))
        end do
#endif
        !        call Driver_abortFlash("DENS var exceeding acceptable range")
     end if
#endif

#ifdef ENER_VAR
     ! energy
     if (any(solnData(ENER_VAR,il:iu,jl:ju,kl:ku) .LT. gr_smalle*0.999999999)) then
        write (block_buff(1,1), '(a)') 'min. solnData(ENER_VAR)'
        write (number_to_str, '('//REAL_FORMAT//')') minval(solnData(ENER_VAR,il:iu,jl:ju,kl:ku))
        write (block_buff(1,2), '(a)') trim(adjustl(number_to_str))

        write (block_buff(2,1), '(a)') 'PE'
        write (number_to_str, '(i7)') myPe
        write (block_buff(2,2), '(a)') trim(adjustl(number_to_str))

        write (block_buff(3,1), '(a)') 'block'
        write (number_to_str, '(i7)') blockID
        write (block_buff(3,2), '(a)') trim(adjustl(number_to_str))

        write (block_buff(4,1), '(a)') 'type'
        write (number_to_str, '(i7)') LEAF
        write (block_buff(4,2), '(a)') trim(adjustl(number_to_str))

        call Logfile_stamp(block_buff, 4, 2, 'WARNING '//info)
        print 111, 'WARNING ',info,':', ((block_buff(i,j),j=1,2),i=1,4)
#ifdef DEBUG_CONSCONV
        do j=blkLimitsGC(HIGH,JAXIS), blkLimitsGC(LOW,JAXIS),-1
           print 112, j, (solnData(ENER_VAR,i,j,blkLimits(LOW,KAXIS)), i=blkLimitsGC(LOW,IAXIS),blkLimitsGC(HIGH,IAXIS)) 
        end do
#endif
        !call Driver_abortFlash("ENER var exceeding acceptable range")
     end if
#endif
#ifdef EINT_VAR
     if (any(solnData(EINT_VAR,il:iu,jl:ju,kl:ku) .LT. gr_smalle*0.999999999)) then
        write (block_buff(1,1), '(a)') 'min. solnData(EINT_VAR)'
        write (number_to_str, '('//REAL_FORMAT//')') minval(solnData(EINT_VAR,il:iu,jl:ju,kl:ku))
        write (block_buff(1,2), '(a)') trim(adjustl(number_to_str))

        write (block_buff(2,1), '(a)') 'PE'
        write (number_to_str, '(i7)') myPe
        write (block_buff(2,2), '(a)') trim(adjustl(number_to_str))

        write (block_buff(3,1), '(a)') 'block'
        write (number_to_str, '(i7)') blockID
        write (block_buff(3,2), '(a)') trim(adjustl(number_to_str))

        write (block_buff(4,1), '(a)') 'type'
        write (number_to_str, '(i7)') LEAF
        write (block_buff(4,2), '(a)') trim(adjustl(number_to_str))

        call Logfile_stamp(block_buff, 4, 2, 'WARNING '//info)
        print 111, 'WARNING ',info,':', ((block_buff(i,j),j=1,2),i=1,4)
#ifdef DEBUG_CONSCONV
        do j=blkLimitsGC(HIGH,JAXIS), blkLimitsGC(LOW,JAXIS),-1
           print 112, j, (solnData(EINT_VAR,i,j,blkLimits(LOW,KAXIS)), i=blkLimitsGC(LOW,IAXIS),blkLimitsGC(HIGH,IAXIS)) 
        end do
#endif
        !          call Driver_abortFlash("EINT var exceeding acceptable range")
     end if
#endif

     call Grid_releaseBlkPtr(blockID, solnData)

  end do

  return 
end subroutine gr_sanitizeDataAfterInterp
