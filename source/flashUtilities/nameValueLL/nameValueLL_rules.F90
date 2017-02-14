!!****ih* source/flashUtilities/nameValueLL/nameValueLL_rules
!!
!!  NAME
!!    nameValueLL_rules
!!
!!  SYNOPSIS
!!
!!    nameValueLL_rules(               type(context_type) :: context,
!!                           character(len=*), intent(in) :: name,
!!                           integer(in)                  :: numValues,
!!                           integer, dimension(numValues):: minValues,
!!                           integer, dimension(numValues):: maxValues)
!!
!!  DESCRIPTION
!!    sets the set of valid values of "name" Runtimeparameter to validValues
!!
!!  ARGUMENTS
!!
!!    context--      data structure holding Runtime Parameters details
!!    name--         name of parameter
!!    numValues--    number of valid values
!!    minValues--    array of given size, provides minimum of valid values
!!    maxValues--    array of given size, provides maximum of valid values
!!
!!***

subroutine nameValueLL_rulesInt ( context, name, numValues, minValues, maxValues)

  use nameValueLL_data !, ONLY: context_type, nameValueLL_findInt, &
 !     &  name_invalid, name_real, name_int, name_str, name_log, &
 !     &  real_list_type, int_list_type, str_list_type, log_list_type
  use Driver_interface, ONLY : Driver_abortFlash

  implicit none

   type (context_type), intent(inout) :: context
   character(len=*), intent(in) :: name
   integer, intent(in) :: numValues
   integer,dimension(numValues),intent(in):: minValues,maxValues
   type (int_list_type), pointer :: node
   integer :: istat

   call nameValueLL_findInt(context,name,node)

   if (.not. associated(node)) then
      call Driver_abortFlash("rules: add parameter and then the rules")
   endif

   node%numValues = 0
   nullify(node%maxValues)
   nullify(node%minValues)
   allocate(node%minValues(numValues),stat=istat)
   if (istat /= 0) then
      call Driver_abortFlash("rules: unable to allocate")
   endif
   allocate(node%maxValues(numValues),stat=istat)
   if (istat /= 0) then
      call Driver_abortFlash("rules: unable to allocate")
   endif
   node%minValues = minValues
   node%maxValues = maxValues
   node%numValues = numValues

end subroutine nameValueLL_rulesInt
   
!!****ih* source/flashUtilities/nameValueLL/nameValueLL_rulesReal
!!
!!  NAME
!!    nameValueLL_rulesReal
!!
!!  SYNOPSIS
!!
!!    nameValueLL_rulesReal ( type(context_type) :: context,
!!                           character(len=*), intent(in) :: name,
!!                           integer(in) :: numValues,
!!                           real, dimension(numValues):: minValues,
!!                           real, dimension(numValues):: maxValues)
!!
!!  DESCRIPTION
!!    sets the set of validvalues of "name" Runtimeparameter to validValues
!!
!!  Arguments
!!    context: data structure holding Runtime Parameters details
!!    name:  name of parameter
!!    numValues: # of valid values
!!    minValues,maxValues: array of values of given size
!!
!!***

subroutine nameValueLL_rulesReal ( context, name, numValues, minValues, maxValues)

  use nameValueLL_data !, ONLY: context_type, nameValueLL_findReal, &
 !     &  name_invalid, name_real, name_int, name_str, name_log, &
 !     &  real_list_type, int_list_type, str_list_type, log_list_type
  use Driver_interface, ONLY : Driver_abortFlash

  implicit none

   type (context_type), intent(inout) :: context
   character(len=*), intent(in) :: name
   integer,intent(in) :: numValues
   real,dimension(numValues),intent(in):: minValues,maxValues
   type (real_list_type), pointer :: node
   integer :: istat

   call nameValueLL_findReal(context,name,node)

   if (.not. associated(node)) then
      call Driver_abortFlash("rules: add parameter and then the rules")
   endif

   node%numValues = 0
   nullify(node%minValues)
   nullify(node%maxValues)
   allocate(node%minValues(numValues),stat=istat)
   if (istat /= 0) then
      call Driver_abortFlash("rules: unable to allocate")
   endif
   allocate(node%maxValues(numValues),stat=istat)
   if (istat /= 0) then
      call Driver_abortFlash("rules: unable to allocate")
   endif
   node%minValues = minValues
   node%maxValues = maxValues
   node%numValues = numValues

end subroutine nameValueLL_rulesReal
   
!!****ih* source/flashUtilities/nameValueLL/nameValueLL_rulesStr
!!
!!  NAME
!!    nameValueLL_rulesStr
!!
!!  SYNOPSIS
!!
!!    nameValueLL_rulesStr ( type(context_type) :: context,
!!                           character(len=*), intent(in) :: name,
!!                           integer(in) :: numValues,
!!                           character(len=MAX_STRING_LENGTH),dimension(numValues):: validValues)
!!
!!  DESCRIPTION
!!    sets the set of validvalues of "name" Runtimeparameter to validValues
!!
!!  Arguments
!!    context: data structure holding Runtime Parameters details
!!    name:  name of parameter
!!    numValues: # of valid values
!!    validValues: array of valid values of given size
!!
!!***

subroutine nameValueLL_rulesStr ( context, name, numValues, validValues)

  use nameValueLL_data !, ONLY: context_type, nameValueLL_findStr, &
 !     &  name_invalid, name_real, name_int, name_str, name_log, &
 !     &  real_list_type, int_list_type, str_list_type, log_list_type
  use Driver_interface, ONLY : Driver_abortFlash

  implicit none

   type (context_type), intent(inout) :: context
   character(len=*), intent(in) :: name
   integer,intent(in) :: numValues
   character(len=*),dimension(numValues),intent(in):: validValues
   type (str_list_type), pointer :: node
   integer :: istat

   call nameValueLL_findStr(context,name,node)

   if (.not. associated(node)) then
      call Driver_abortFlash("rules: add parameter and then the rules")
   endif

   node%numValues = 0
   nullify(node%validValues)
   allocate(node%validValues(numValues),stat=istat)
   if (istat /= 0) then
      call Driver_abortFlash("rules: unable to allocate")
   endif
   node%validValues = validValues
   node%numValues = numValues
   do istat = 1 , numValues
      call makeLowercase(node%validValues(istat))
   enddo

end subroutine nameValueLL_rulesStr

