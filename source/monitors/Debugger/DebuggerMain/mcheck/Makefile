# The +mcheck setup shortcut adds MCHECK macro to CDEFINES.
#
# The Makefile fragment below modifies LIB if CDEFINES contains
# MCHECK.  Manipulation of LIB in this way is extremely ugly, but it
# ensures that the malloc replacement library is last in the link
# line.  This cannot be guaranteed if using the library capabililty of
# FLASH.

ifneq (,$(findstring MCHECK,$(CDEFINES)))
    LIB += -lmcheck
endif
