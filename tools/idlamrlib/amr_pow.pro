function amr_pow, amr1, ex

; Copyright Mark Krumholz and Robert Fisher (2001)
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


; this function does a cell-wise exponentation of amr1, raising it to
; the power pow, and returning a new amr object

; Create the new amr object to hold the result. Null out the pointers
; in it.
amr_new = amr1
for n=0, amr_new.maxlevel do begin
	for m=0, amr_new.levels[n].nfab do $
		amr_new.levels[n].fabptr = ptr_new()
endfor

; change component name field in new object
amr_new.componentname = '(' + amr1.componentname + '^' + $
	strtrim(string(ex),2) + ')'

; go through levels
for n=0, amr1.maxlevel do begin

	; create a new fab pointer for amr_new
	amr_new.levels[n].fabptr = ptr_new(*amr1.levels[n].fabptr)

	; for safety, null out all the data pointers within this fab pointer
	for m=0, amr_new.levels[n].nfab-1 do $
		(*amr_new.levels[n].fabptr)[m].dataptr = ptr_new()

	; now go through fabs
	for m=0, amr1.levels[n].nfab-1 do begin

		; form the power, and set the new data pointer to point to it
		datanew = (*(*amr1.levels[n].fabptr)[m].dataptr)^ex
		(*amr_new.levels[n].fabptr)[m].dataptr = ptr_new(datanew)

	endfor
endfor

; return the new amr object
return, amr_new
end

