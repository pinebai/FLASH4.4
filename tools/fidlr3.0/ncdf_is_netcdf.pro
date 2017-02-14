FUNCTION ncdf_is_netcdf, FILENAME

;  A function to check if a filename is truly in netCDF version
;  A function similar to H5F_IS_HDF5
;  with great thanks to
;  http://astrog.physics.wisc.edu/~craigm/idl/archive/msg00581.html
;  cleaned up by Lynn Reid 14/09/2005

;- Set return values

false = 0B
true = 1B

;- Establish error handler
catch, error_status
IF (error_status ne 0) THEN BEGIN 
  catch, /cancel
  return, false
ENDIF 

;- Try opening the file
cdfid = NCDF_OPEN( filename,/NOWRITE )

;- If we get this far, open must have worked
NCDF_CLOSE, cdfid
catch, /cancel
return, true

END
