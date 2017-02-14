!!****if* source/physics/Hydro/HydroMain/split/PPM/PPMKernel/hy_nomonot
!!
!! NAME
!! 
!!  hy_nomonot
!!
!! SYNOPSIS
!!
!!  call hy_nomonot(integer(IN) :: numIntCells, 
!!             integer(IN) :: numCells, 
!!             real(INOUT), dimension(numCells) :: al,
!!             real(IN),    dimension(numCells) :: a,
!!             real(INOUT), dimension(numCells) :: ar,
!!             real(OUT),   dimension(numCells) :: da,
!!             real(OUT),   dimension(numCells) :: a6)
!!
!! 
!! DESCRIPTION
!!
!!  Apply no monotonicity constraint to interpolation parabola -- do not constrain
!!  the parabolic distribution of each variable so that all points in the
!!  parabola fall between the zone interface values.
!!
!!
!! ARGUMENTS
!!   numIntCells - number if interior cells
!!   numCells    - maximum of the number of cells along any dimension
!!   al          - array containing left edge values
!!   a           - array containing center values
!!   ar          - array containing right edge values
!!   da          - derivative values
!!   a6          - 
!!
!! SEE ALSO
!!
!!  monot
!!***

subroutine hy_nomonot(numIntCells, numCells,al,a,ar,da,a6)

  use Hydro_data, ONLY: hy_iplm

  implicit none
  integer, intent(IN) :: numIntCells, numCells
  real, intent(INOUT), DIMENSION(numCells) :: al, ar
  real, intent(IN),    DIMENSION(numCells) :: a
  real, intent(OUT),   DIMENSION(numCells) :: da, a6
  
  real :: disval
  integer :: i,numIntCells5
  real s1, s2, s3
  
  numIntCells5 = numIntCells + 5
  
  do  i = 4, numIntCells5
     
     da(i) = (ar(i) - al(i))
     
     s1 = (ar(i) - a(i)) * (al(i) - a(i))
     
!!$     if (s1 >= 0.e0) then
!!$        al(i) = a(i)
!!$        ar(i) = a(i)
!!$     endif
     

! Eq. 1.10 of Colella & Woodward
    
     disval = (ar(i) - a(i)) * (al(i) - a(i))
     
!!$     if (disval /= 0.e0) then
        s2 = 3.e00 * a(i) - 2.e00 * ar(i)
        s3 = 3.e00 * a(i) - 2.e00 * al(i)
!!$     else
!!$        s2 = al(i)
!!$        s3 = ar(i)
!!$     endif
       

     da(i) = ar(i) - al(i)


! if we are doing the Piecewise Linear Method, enforce it here

     if (hy_iplm == 0) then
        a6(i) = 6.e00 * a(i) - 3.e00 * (al(i) + ar(i))
     else
        a6(i) = 0.e00
     endif
     
  end do
  
  return
end subroutine hy_nomonot
  



