function determine_flash_version, filename
;
; determine whether a file is flash2 or flash 3
;  essentially the same as determine_file_version, but
; RETURNS:  -1 if file does not exist or is an unknown format
;            2 if FLASH2
;            3 if FLASH3
;

; first determine if the file exists!
file_exists = file_test(filename)
IF (NOT file_exists) THEN BEGIN
  print,'ERROR: file ',filename,' does not exist.'
  return, -1
ENDIF 

; this does the real work
file_version = determine_file_version(filename)

; initialize to unknown
flash_version = -1
IF (file_version EQ 7) THEN flash_version = 2
IF (file_version GE 8) THEN flash_version = 3

return, flash_version


END ; determine_flash_version.pro 

