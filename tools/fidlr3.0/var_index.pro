function var_index, string

; translation function to take the variable name into an index for
; setting contour limits, etc.

; this is a general version, for use in the main plot routines
;
; fidlr uses two arrays to hold the variable names, varnames contains
; all the variable names in the data file, plus any derived variables
; that can be made from them.  unk_names contains the names of the 
; variables in the file.  If varnames is not defined, used unk_names

; common in the variable strings
common variables, varnames, native_num
common vars, unk, unk_names

undefined = 0

if n_elements(varnames) EQ 0 AND n_elements(unk_names) GT 0 then $
  varnames = unk_names

if n_elements(varnames) GT 0 then begin
  if varnames[0] EQ '----' then undefined = 1
endif else begin
  undefined = 1
endelse


if undefined then begin

  ; if varnames is not yet defined, then set the index to the max of the
  ; ctr_lim array, so we don't go out of bounds in xflash_defaults
  index = 99

endif else begin
  index = (where(varnames EQ string))[0]
endelse

return, index

end
