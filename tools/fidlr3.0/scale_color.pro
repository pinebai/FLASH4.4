;------------------------------------------------------------------------------
; scale_color.pro -- MZ 3-10-00
;
; take the hydro variable in block format and create a single 1 byte
; array at the resolution of the finest block for plotting purposes.  
; Optionally allow a resolution less than the finest to create a
; smaller array (use the sample keyword)
;
;
; arguments:
;
;     var --  the variable to be scaled, var[maxblocks,nxb,nyb]
;
; optional arguments
;
;     /log    -- scale to the log of var
;
;     varmin  \  set the minimum and maximum limits of var to scale
;     varmax  /  between
;
;     colormap_min \  The limits of the colormap to use -- default is
;     colormap_max /  0 to 255
;------------------------------------------------------------------------------

function scale_color, var, VARMAX = varmax, VARMIN = varmin, $
                      LOG = log, $
                      COLORMAP_MIN = colormap_min, $
                      COLORMAP_MAX = colormap_max

temp_arr = var

if n_elements(log) EQ 0 then log = 0
if n_elements(varmax) EQ 0 then varmax = max(temp_arr)
if n_elements(varmin) EQ 0 then varmin = min(temp_arr)

if n_elements(colormap_min) EQ 0 then colormap_min = 0
if n_elements(colormap_max) EQ 0 then colormap_max = 255

if log EQ 1 then begin
    varmax   = alog10(varmax)
    varmin   = alog10(varmin)
endif


; mark those regions which are holes -- -1.e30, and make these black
iholes = where(temp_arr EQ -1.e30)

if colormap_min GT colormap_max then begin
    temp_arr = colormap_min - bytscl(temporary(temp_arr), $
                                       min = varmin, max = varmax, $
                                       top = colormap_min - colormap_max)
endif else begin
    temp_arr = colormap_min + bytscl(temporary(temp_arr), $
                                       min = varmin, max = varmax, $
                                       top = colormap_max - colormap_min)
endelse

; mark the holes black
if min(iholes GT 0) then temp_arr[iholes] = 0
    
return, temp_arr

end


