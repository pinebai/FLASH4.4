"""Small utility functions for doing useful stuff for flash."""
import os

def strip_comments(line, comment_char='#', quote_char='"'):
    """Returns the string after removing comments."""
    cut = 0
    while cut < len(line):
       pos = line[cut:].find(comment_char)

       # no real comment in line
       if pos < 0: 
          return line  

       # do we have even number quotes before it? If so this is real comment
       if line[:cut + pos].count(quote_char) % 2 == 0:
          return line[:cut + pos]

       cut = cut + pos + 1
    return line


#########################
### message functions ###
#########################

USE_COLOR = (os.name is 'posix')

def message(s):
    """Formats a message for printing.  If on a posix system the message will be in color."""
    head = "\033[1;32m" if USE_COLOR else "*** MESSAGE ***: "
    tail = "\033[0m" if USE_COLOR else ""
    msg = head + s + tail
    return msg


def warning(s):
    """Formats a message for warning.  If on a posix system the message will be in color."""
    head = "\033[1;33m" if USE_COLOR else "*** WARNING ***: "
    tail = "\033[0m" if USE_COLOR else ""
    msg = head + s + tail
    return msg


def failure(s):
    """Formats a fail message for printing.  If on a posix system the message will be in color."""
    head = "\033[1;31m" if USE_COLOR else "*** FAILURE ***: "
    tail = "\033[0m" if USE_COLOR else ""
    msg = head + s + tail
    return msg
