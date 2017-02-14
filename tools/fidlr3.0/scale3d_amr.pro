;------------------------------------------------------------------------------
; scale3d_amr.pro -- MZ 1-18-99
;
; take the hydro variable in block format and create a single 1 byte
; array at the resolution of the finest block for plotting purposes.  
; Optionally allow a resolution less than the finest to create a
; smaller array (use the sample keyword)
;
; arguments:
;
;     var --  the variable to be scaled, var[maxblocks,nxb,nyb]
;
;
;     colormap_min \  The limits of the colormap to use -- default is
;     colormap_max /  0 to 255
;------------------------------------------------------------------------------

function scale3d_amr, var, VARMIN = varmin, VARMAX = varmax

if n_elements(varmax) EQ 0 then varmax = max(var)
if n_elements(varmin) EQ 0 then varmin = min(var)

; mark those regions which are holes -- -1.e30, and make these black
iholes = where(var EQ -1.e30)

s = size(var)
;temp_arr = var

temp_scale = bytarr(s[1],s[2],s[3])

temp_scale = bytscl(var, min = varmin, max = varmax)

; mark the holes black
if min(iholes GT 0) then temp_scale[iholes] = 0
    
return, temp_scale

end


