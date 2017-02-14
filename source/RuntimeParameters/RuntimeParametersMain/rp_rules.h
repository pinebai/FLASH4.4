
interface rp_rules

  subroutine rp_rulesInt (name, numValues, minValues, maxValues)
    character(len=*), intent(in)          :: name
    integer, intent(in)                   :: numValues
    integer, dimension(numValues),intent(in) :: minValues, maxValues
  end subroutine rp_rulesInt

  subroutine rp_rulesReal (name, numValues, minValues, maxValues)
    character(len=*), intent(in)          :: name
    integer, intent(in)                   :: numValues
    real, dimension(numValues),intent(in) :: minValues, maxValues
  end subroutine rp_rulesReal

  subroutine rp_rulesStr (name, numValues, validValues)
    character(len=*), intent(in)          :: name
    integer, intent(in)                   :: numValues
    character(len=*), dimension(numValues), intent(in)          :: validValues
  end subroutine rp_rulesStr

end interface

