Module ut_randomInterface
  interface
     subroutine ut_randomSeed( ut_size, ut_put, ut_get)
       integer, optional, dimension(:), intent(OUT) :: ut_get
       integer, optional, dimension(:), intent(IN) :: ut_put
       integer, optional, intent(OUT) :: ut_size
     end subroutine ut_randomSeed
  end interface

  interface
     subroutine ut_randomNumber(x)
       real, intent(OUT) :: x
     end subroutine ut_randomNumber
  end interface

  interface
     subroutine ut_randomNumberArray(x)
       real, intent(OUT) :: x(:)
     end subroutine ut_randomNumberArray
  end interface

end Module ut_randomInterface
