#!/bin/bash
# The value 7FF7FFFF, when stored in a REAL(4) variable, has a
# single-precision NaNQ value. If the same number is stored twice in a
# REAL(8) variable (7FF7FFFF7FF7FFFF), it has a double-precision NaNS
# value.

# Signal handlers
# ---------------
# xl__ieee
# Produces a traceback and an explanation of the signal and continues
# execution by supplying the default IEEE result for the failed
# computation. This handler allows the program to produce the same
# results as if exception detection was not turned on.
#
# xl__trcedump
# Produces a traceback and a core file and stops the program.

SIGNALING_NAN='7ff7ffff'
QUIET_NAN='ff'

#NAN=${SIGNALING_NAN}
NAN=${QUIET_NAN}

xlf90_r test_fpe_events.F90 -c -g -qrealsize=8 \
    -qinitauto=${NAN} -qinitalloc=${NAN} \
    -qflttrap=enable:invalid:overflow:zerodivide \
    -qsigtrap=xl__ieee
xlf90_r test_fpe_events.o -o test_fpe_events
