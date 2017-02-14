!!****if* source/flashUtilities/nameValueLL/findWords
!!
!! NAME
!!  findWords
!!
!! SYNOPSIS
!!
!!  findWords(character(len=*)(INOUT) :: string,
!!            integer(INOUT)          :: first(n),
!!            integer(INOUT)          :: last(n),
!!            integer(IN)             :: n,
!!            integer(INOUT)          :: nwords)
!!
!! DESCRIPTION
!!
!!    Given an input string, return the character positions of the
!!    beginnings and endings of all words in the string.  Also
!!    return the number of words found.
!!
!!    What constitutes a word:
!!
!!       <delimiter>characters<delimiter>
!!       symbol (single character)
!!       <delimiter>characters<end-of-line>
!!       <beg-of-line>characters<delimiter>
!!
!!        Symbols also function as delimiters.
!!
!!
!! ARGUMENTS
!!
!!  string  -- input string
!!  first(n)-- array of character positions of the beginning of all words in string
!!  last(n) -- array of character positions of the end of all words in string
!!  n       -- size of arrays first and last
!!  nwords  -- number of words found in the string
!!
!!***
      
subroutine findWords (string, first, last, n, nwords)
  
  implicit none

#include "constants.h"
  
  character(len=*), intent(inout)   :: string
  integer, intent(in)               :: n 
  integer, intent(inout)            :: nwords, first(n), last(n)
  
  integer, parameter :: ndel = 2, nsym = 1
  !The deliminator characters must be a space and a TAB.
  character    :: delimiters(ndel) = (/' ', TAB_CHAR/)
  character    :: symbols(nsym)    = (/'='/)
  integer              :: i, j, k, strlen
  logical              :: is_delimiter, is_symbol

  
  
  
  strlen = len(string)
  nwords = 0
  k      = 1
  
  do i = 1, strlen
     
     is_delimiter = .false.
     do j = 1, ndel
        if (string(i:i) == delimiters(j)) is_delimiter = .true.
     enddo
     is_symbol = .false.
     do j = 1, nsym
        if (string(i:i) == symbols(j)) is_symbol = .true.
     enddo
     
     if (is_delimiter .or. is_symbol) then
        if (i > k) then
           nwords = nwords + 1
           first(nwords) = k
           last(nwords)  = i - 1
        endif
        k = i + 1
     endif
     if (is_symbol) then
        nwords = nwords + 1
        first(nwords) = i
        last(nwords)  = i
     endif
     
  enddo
  
  !               If we ended without a delimiter or symbol, record the last
  !               word.
  
  if (.not. (is_delimiter .or. is_symbol)) then
     nwords = nwords + 1
     first(nwords) = k
     last(nwords)  = strlen
  endif
  
  
  return
end subroutine findWords


