function amr_log, amr1, base

; Copyright Robert Fisher (2001)
;
; This program is free software; you can redistribute it and/or modify
; it under the terms of the GNU General Public License as published by
; the Free Software Foundation; either version 2 of the License, or
; (at your option) any later version.
;
; This program is distributed in the hope that it will be useful,
; but WITHOUT ANY WARRANTY; without even the implied warranty of
; MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
; GNU General Public License for more details.
;
; You should have received a copy of the GNU General Public License
; along with this program; if not, write to the Free Software
; Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

; this function does a cell-wise log of amr1, taken to the given base
; (default is 10)

; Set base
if (n_elements(base) eq 0) then base = 10.0
if (base le 0) then begin
	message, 'Error: base must be positive!'
	return, -1
endif

; Create the new amr object to hold the result. Null out the pointers
; in it.
amr_new = amr1
for n=0, amr_new.maxlevel do begin
        for m=0, amr_new.levels[n].nfab do $
                amr_new.levels[n].fabptr = ptr_new()
endfor

; change component name field in new object
amr_new.componentname = '(' + amr1.componentname + 'log)' 

; go through levels
for n=0, amr1.maxlevel do begin

	; create a new fab pointer for amr_new
	amr_new.levels[n].fabptr = ptr_new(*amr1.levels[n].fabptr)

	; for safety, null out all the data pointers within this fab pointer
	for m=0, amr_new.levels[n].nfab-1 do $
		(*amr_new.levels[n].fabptr)[m].dataptr = ptr_new()

	; now go through fabs
	for m=0, amr1.levels[n].nfab-1 do begin

	; check for non-positive values
	if total(where((*(*amr1.levels[n].fabptr)[m].dataptr) le 0.0)) ne -1 $
		then message, 'Warning: non-positive cell value detected' + $
		' in amr_log, using absolute value', /continue

	; form the log, and set the new data pointer to point to it
	datanew = alog (abs (*(*amr1.levels[n].fabptr)[m].dataptr) ) / $
		alog(base) 
	(*amr_new.levels[n].fabptr)[m].dataptr = ptr_new(datanew)

        endfor
endfor

; return the new amr object
return, amr_new
end

