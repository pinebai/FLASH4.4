function coord_to_val, coord, amr, interp=interp

; Derived from coord_to_fab (C) Mark Krumholz

; This program returns the value of an amr object at the coordinate
; specified
; Setting Interp true makes it use [bi|tri]linear interpolation
; The default behavior is to sample the nearest cell center 
;reslevels = 0

; go through amr structure
res = !values.f_nan
for n=amr.maxlevel, 0,-1 do begin
    for m=0, amr.levels[n].nfab-1 do begin
        
        intersect = 1
        
        if (coord[0] lt (*amr.levels[n].fabptr)[m].xlo[0]) or $
          (coord[0] gt (*amr.levels[n].fabptr)[m].xhi[0]) then $
          intersect = 0
        if amr.ndim ge 2 then begin
            if (coord[1] lt (*amr.levels[n].fabptr)[m].xlo[1]) or $
              (coord[1] gt (*amr.levels[n].fabptr)[m].xhi[1]) then $
              intersect = 0
        endif
        if amr.ndim ge 3 then begin
            if (coord[2] lt (*amr.levels[n].fabptr)[m].xlo[2]) or $
              (coord[2] gt (*amr.levels[n].fabptr)[m].xhi[2]) then $
              intersect = 0
        endif
        
        if not intersect then continue else begin
            ; It intersects... find the value
            if keyword_set(interp) then begin
                print,'Interp keyword not yet implemented for COORD_TO_VAL'
            endif else begin
                offset = coord - (*amr.levels[n].fabptr)[m].xlo
                fabidx = round((offset - 0.5*amr.gridspacing[*,n]) / $
                               amr.gridspacing[*,n])
                fabidx = fabidx * (fabidx ne -1)
                fabidxmax = (*amr.levels[n].fabptr)[m].idxhi - $
                  (*amr.levels[n].fabptr)[m].idxlo
                fabidx = fabidx * (fabidx ne fabidxmax) + $
                  fabidxmax * (fabidx eq fabidxmax)
                if amr.ndim eq 1 then res = (*((*amr.levels[n].fabptr)[m]).dataptr)[fabidx[0]]
                if amr.ndim eq 2 then res = (*((*amr.levels[n].fabptr)[m]).dataptr)[fabidx[0], fabidx[1]]
                if amr.ndim eq 3 then res = (*((*amr.levels[n].fabptr)[m]).dataptr)[fabidx[0], fabidx[1], fabidx[2]]
            endelse
            m = amr.levels[n].nfab
            n = -1
        endelse
    endfor
endfor

return,res

;if reslevels ne 0 then return, res else return, -1

end
