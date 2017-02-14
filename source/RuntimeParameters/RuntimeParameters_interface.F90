!!****h* source/RuntimeParameters/RuntimeParameters_interface
!!
!! This is the header file for the RuntimeParameters module
!! that defines its public interfaces.
!!***
Module RuntimeParameters_interface

  implicit none  

# include "Flash.h"
# include "constants.h"

  interface RuntimeParameters_add
    subroutine RuntimeParameters_addReal (name, value, rwState)
      character(len=*), intent(in)          :: name
      real, intent(in)                      :: value
      integer,OPTIONAL,intent(in)           :: rwState
    end subroutine RuntimeParameters_addReal

    subroutine RuntimeParameters_addInt (name, value, rwState)
      character(len=*), intent(in)          :: name
      integer, intent(in)                   :: value
      integer,OPTIONAL,intent(in)           :: rwState
    end subroutine RuntimeParameters_addInt

    subroutine RuntimeParameters_addStr (name, value, rwState)
      character(len=*), intent(in)          :: name
      character(len=*), intent(in)          :: value
      integer,OPTIONAL,intent(in)           :: rwState
    end subroutine RuntimeParameters_addStr

    subroutine RuntimeParameters_addLog (name, value, rwState)
      character(len=*), intent(in)          :: name
      logical, intent(in)                   :: value
      integer,OPTIONAL,intent(in)           :: rwState
    end subroutine RuntimeParameters_addLog
  end interface

  interface
    subroutine RuntimeParameters_bcast()
    end subroutine RuntimeParameters_bcast
  end interface

  interface
    subroutine RuntimeParameters_finalize()
    end subroutine RuntimeParameters_finalize
  end interface

  interface RuntimeParameters_get
    subroutine RuntimeParameters_getReal (name, value)
      character(len=*), intent(in)          :: name
      real, intent(out)                     :: value
    end subroutine RuntimeParameters_getReal

    subroutine RuntimeParameters_getInt (name, value)
      character(len=*), intent(in)          :: name
      integer, intent(out)                  :: value
    end subroutine RuntimeParameters_getInt

    subroutine RuntimeParameters_getStr (name, value)
      character(len=*),intent(in)           :: name
      character(len=*),intent(out)          :: value
    end subroutine RuntimeParameters_getStr

    subroutine RuntimeParameters_getLog (name, value)
      character(len=*),intent(in)           :: name
      logical,intent(out)                   :: value
    end subroutine RuntimeParameters_getLog
  end interface

  interface RuntimeParameters_getAll
    subroutine RuntimeParameters_getAllReal (num, names, values, changed)
      integer, intent(in)                                   :: num
      character(len=MAX_STRING_LENGTH), intent(inout)       :: names(num)
      real, intent(inout)                                   :: values(num)
      logical, intent(inout)                                :: changed(num)
    end subroutine RuntimeParameters_getAllReal

    subroutine RuntimeParameters_getAllInt (num, names, values, changed)
      integer, intent(in)                                   :: num
      character(len=MAX_STRING_LENGTH), intent(inout)       :: names(num)
      integer, intent(inout)                                :: values(num)
      logical, intent(inout)                                :: changed(num)
    end subroutine RuntimeParameters_getAllInt

    subroutine RuntimeParameters_getAllStr (num, names, values, changed)
      integer, intent(inout)                                :: num
      character(len=MAX_STRING_LENGTH), intent(inout)       :: names(num), values(num)
      logical, intent(inout)                                :: changed(num)
    end subroutine RuntimeParameters_getAllStr

    subroutine RuntimeParameters_getAllLog (num, names, values, changed)
      integer, intent(in)                                     :: num
      character(len=MAX_STRING_LENGTH), intent(inout)         :: names(num)
      logical, intent(inout)                                  :: values(num)
      logical, intent(inout)                                  :: changed(num)
    end subroutine RuntimeParameters_getAllLog
  end interface

  ! Even though the API implies there is an overloaded routine call
  !  RuntimeParameters_getNum, there is NOT.  The routines have
  !  to be called directly 
  interface
    subroutine RuntimeParameters_getNumReal (nparms)
      integer, intent(out)                :: nparms
    end subroutine RuntimeParameters_getNumReal
  end interface

  interface
    subroutine RuntimeParameters_getNumInt (nparms)
      integer, intent(out)               :: nparms
    end subroutine RuntimeParameters_getNumInt
  end interface

  interface
    subroutine RuntimeParameters_getNumStr (nparms)
      integer, intent(out)               :: nparms
    end subroutine RuntimeParameters_getNumStr
  end interface

  interface
    subroutine RuntimeParameters_getNumLog (nparms)
      integer, intent(out)               :: nparms
    end subroutine RuntimeParameters_getNumLog
  end interface

  interface
    subroutine RuntimeParameters_getNumIgn (numIgnoredParams)
      integer, intent(out)               :: numIgnoredParams
    end subroutine RuntimeParameters_getNumIgn
  end interface

  interface RuntimeParameters_getPrev
    subroutine RuntimeParameters_getPrevReal (name, value)
      character(len=*), intent(in)          :: name
      real, intent(out)                     :: value
    end subroutine RuntimeParameters_getPrevReal

    subroutine RuntimeParameters_getPrevInt (name, value)
      character(len=*), intent(in)          :: name
      integer, intent(out)                  :: value
    end subroutine RuntimeParameters_getPrevInt

    subroutine RuntimeParameters_getPrevStr (name, value)
      character(len=*),intent(in)           :: name
      character(len=*),intent(out)          :: value
    end subroutine RuntimeParameters_getPrevStr

    subroutine RuntimeParameters_getPrevLog (name, value)
      character(len=*),intent(in)           :: name
      logical,intent(out)                   :: value
    end subroutine RuntimeParameters_getPrevLog
  end interface

  interface
    subroutine RuntimeParameters_init(restart)
      logical, intent(out)                   :: restart
    end subroutine RuntimeParameters_init
  end interface

  interface
    subroutine RuntimeParameters_mapStrToInt (inputString, constKey)
      character(len=MAX_STRING_LENGTH), intent(in) :: inputString
      integer, intent(inout)                :: constKey
    end subroutine RuntimeParameters_mapStrToInt
  end interface

  interface
    subroutine RuntimeParameters_read (parmfile)
      character(len=MAX_STRING_LENGTH), intent(in)    :: parmfile
    end subroutine RuntimeParameters_read
  end interface

  interface RuntimeParameters_set
    subroutine RuntimeParameters_setReal (name, value)
      character(len=*), intent(in)          :: name
      real, intent(in)                      :: value
    end subroutine RuntimeParameters_setReal

    subroutine RuntimeParameters_setInt (name, value)
      character(len=*), intent(in)          :: name
      integer, intent(in)                   :: value
    end subroutine RuntimeParameters_setInt

    subroutine RuntimeParameters_setStr (name, value)
      character(len=*),intent(in)           :: name, value
    end subroutine RuntimeParameters_setStr

    subroutine RuntimeParameters_setLog (name, value)
      character(len=*),intent(in)           :: name
      logical,intent(in)                    :: value
    end subroutine RuntimeParameters_setLog
  end interface

  interface RuntimeParameters_setPrev
    subroutine RuntimeParameters_setPrevReal (name, value)
      character(len=*), intent(in)          :: name
      real, intent(in)                      :: value
    end subroutine RuntimeParameters_setPrevReal
  
    subroutine RuntimeParameters_setPrevInt (name, value)
      character(len=*), intent(in)          :: name
      integer, intent(in)                   :: value
    end subroutine RuntimeParameters_setPrevInt

    subroutine RuntimeParameters_setPrevStr (name, value)
      character(len=*),intent(in)           :: name, value
    end subroutine RuntimeParameters_setPrevStr

    subroutine RuntimeParameters_setPrevLog (name, value)
      character(len=*),intent(in)           :: name
      logical,intent(in)                    :: value
    end subroutine RuntimeParameters_setPrevLog
  end interface

  interface
    subroutine RuntimeParameters_stampIgnored()
    end subroutine RuntimeParameters_stampIgnored
  end interface

end Module RuntimeParameters_interface
