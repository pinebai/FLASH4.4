pro get_lrefine_max_min, filename, ndim, MIN=lrefine_min, MAX=lrefine_max

; given a filename, get the max and min refinement levels

if (n_elements(filename) EQ 0) then begin
    print, 'ERROR: no filename specified to get_lrefine_max_min'
    return
endif

itype = determine_file_type(filename)

if(itype EQ 1) then begin
;------------------------------------------------------------------------------
; open up the file for read now and read in the refine levels
;------------------------------------------------------------------------------
    
    file_identifier = H5F_OPEN(filename)
    dataset = H5D_OPEN(file_identifier, "refine level")
    lrefine = H5D_READ(dataset)
    H5D_CLOSE, dataset
        
    lrefine_min = min(lrefine)
    lrefine_max = max(lrefine)
    H5F_CLOSE, file_identifier

endif

if(itype EQ 2) then begin

	file_identifier = NCDF_OPEN(filename)
	ncdf_varget, file_identifier, "lrefine", lrefine
	lrefine_min = min(lrefine)
	lrefine_max = max(lrefine)
	NCDF_CLOSE, file_identifier
endif

return

end
