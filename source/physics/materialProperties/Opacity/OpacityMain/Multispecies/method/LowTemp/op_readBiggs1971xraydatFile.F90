!!****if* source/physics/materialProperties/Opacity/OpacityMain/Multispecies/method/LowTemp/op_readBiggs1971xraydatFile
!!
!! NAME
!!
!!  op_readBiggs1971xraydatFile
!!
!! SYNOPSIS
!!
!!  call op_readBiggs1971xraydatFile ()
!!
!! DESCRIPTION
!!
!!  Reads the 1971 F.Biggs and R.Lighthill 'Biggs1971xraydat.txt' file into 1d arrays and produces
!!  array initialization output for use in the opacity unit in FLASH. This routine is just a tool
!!  to produce the appropriate subroutines for FLASH.
!!
!! ARGUMENTS
!!
!!***
subroutine op_readBiggs1971xraydatFile ()

  use Driver_interface, ONLY : Driver_abortFlash

  implicit none

  character (len=80) :: dummyLine

  logical :: fileExists

  integer :: fileUnit, ut_getFreeFileUnit
  integer :: i,j,L,n,nj

  real    :: a1,a2,a3,a4

  integer :: J_beg    (1:101)

  real    :: A_linear (1:3884)
  real    :: EG_lower (1:971)
  real    :: EG_upper (1:971)
  real    :: Z_over_A (1:100)

  real    :: A (4,12,100)
!
!
!   ...Check and open the opacity file.
!
!
  inquire (file = "Biggs1971xraydat.txt" , exist = fileExists)

  if (.not.fileExists) then
       call Driver_abortFlash ('[op_readBiggs1971xraydatFile] ERROR: Biggs1971xraydat.txt file not found')
  end if

  fileUnit = ut_getFreeFileUnit ()
  open (unit = fileUnit , file = "Biggs1971xraydat.txt")
!
!
!   ...Read the info of the xraydat file into the linear arrays.
!
!
  read (fileUnit,'(A80)')  dummyLine
  read (fileUnit,'(A80)')  dummyLine
  read (fileUnit,'(A80)')  dummyLine
  read (fileUnit,'(A80)')  dummyLine
  read (fileUnit,'(A80)')  dummyLine        ! documentation junk that is not needed
  read (fileUnit,'(A80)')  dummyLine
  read (fileUnit,'(A80)')  dummyLine
  read (fileUnit,'(A80)')  dummyLine
  read (fileUnit,'(A80)')  dummyLine
  read (fileUnit,'(A80)')  dummyLine

  read (fileUnit,'(8E10.3)') (A_linear (n), n = 1,3884)          ! the A(i,j,4) expansion coefficients
  read (fileUnit,'(20I4)')   (J_beg    (n), n = 1,100)           ! the 'j' interval beginnings for all 100 atoms 
  read (fileUnit,'(8F10.4)') (EG_upper (n), n = 1,971)           ! the upper radiation energy level for each 'j' interval 
  read (fileUnit,'(11F7.4)') (Z_over_A (n), n = 1,100)           ! the upper radiation energy level for each 'j' interval

  J_beg (101) = 971

  EG_lower (2:971) = EG_upper (1:970)

  L = 1
  do i = 1,100
     nj = J_beg (i+1) - J_beg (i)
     EG_lower (L) = 0.01
     L = L + nj
  end do
!
!
!   ...Write out the EG_lower and EG_upper arrays.
!
!
  L = 0

  do i = 1,100
     nj = J_beg (i+1) - J_beg (i)

     write (*,*)

     do j = 1,nj

        L = L + 1
        write (*,8000) '  op_PEenergyRange (LOW:HIGH,',j,',',i,') = (/',EG_lower (L),',',EG_upper (L),'/)'

     end do
  end do

  8000 format (A29,I2,A1,I3,A6,F10.4,A1,F10.4,A2)
!
!
!   ...Write out the A(4 j,i) array.
!
!
!  L = 1
!
!  do i = 1,100
!     nj = J_beg (i+1) - J_beg (i)
!
!     write (*,*)
!
!     do j = 1,nj
!
!        a1 = A_linear (L)
!        a2 = A_linear (L+1)
!        a3 = A_linear (L+2)
!        a4 = A_linear (L+3)
!
!        write (*,9000) '  op_Aij4 (1:4,',j,',',i,') = (/',a1,',',a2,',',a3,',',a4,'/)'
!
!        L = L + 4
!
!     end do
!  end do
!
!  9000 format (A15,I2,A1,I3,A6,ES10.3,A1,ES10.3,A1,ES10.3,A1,ES10.3,A2)
!
!
!   ...Close the xraydat file.
!
!
  close (fileUnit)
!
!
!   ...Ready! 
!
!
  return
end subroutine op_readBiggs1971xraydatFile
