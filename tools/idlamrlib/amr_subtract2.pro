function amr_subtract2, amr1, amr2

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


; This function does a cell-wise subtraction of two amr objects of
; identical structure and returns a new amr object. If amr2 is a pure
; number rather than an amr object, the routine subtracts amr2 from
; every cell of amr1.

; see whether amr2 is an amr objects or a pure number
if n_tags(amr2) eq 0 then isAMR2 = 0 else isAMR2 = 1

if isAMR2 then begin

	; Start with an error check to make sure the structures are compatible
	; at the highest level.
	if (amr1.ndim ne amr2.ndim) or (amr1.maxlevel ne amr2.maxlevel) or $
		(total(amr1.boxmin ne amr2.boxmin) ne 0) or $
		(total(amr1.boxmax ne amr2.boxmax) ne 0) or $
		(total(amr1.refratio ne amr2.refratio) ne 0) or $
		(total(amr1.idxlo ne amr2.idxlo) ne 0) or $
		(total(amr1.idxhi ne amr2.idxhi) ne 0) or $
		(total(amr1.gridspacing ne amr2.gridspacing) ne 0) or $
		(amr1.coordtype ne amr2.coordtype) then begin
;			print, 'Error: structures are incompatible'
;			return, -1
                  mergestructures_amr, amr1, amr2, outamr1, outamr2
                  amr1 = outamr1
                  amr2 = outamr2
	endif	
endif

; Create the new amr object to hold the result. Null out the pointers
; in it.
amr_new = amr1
for n=0, amr_new.maxlevel do begin
	for m=0, amr_new.levels[n].nfab do $
		amr_new.levels[n].fabptr = ptr_new()
endfor

; change some fields in the new object
if isAMR2 then begin
	if (amr1.name ne amr2.name) then amr_new.name = amr1.name + ' - ' $
		+ amr2.name
	amr_new.componentname = '(' + amr1.componentname + ' - ' $
		+ amr2.componentname + ')'
	if amr1.time ne amr2.time then amr_new.time = -1
	if total(amr1.periodicity ne amr2.periodicity) ne 0 then $
		amr_new.periodicity = -1
	if total(amr1.levelsteps ne amr2.levelsteps) ne 0 then $
		amr_new.levelsteps = -1
endif else begin
	amr_new.componentname = '(' + amr1.componentname + $
		' - ' + strtrim(string(amr2),2) + ')'
endelse

; go through levels
for n=0, amr1.maxlevel do begin

	; do consistency checks on amr1 and amr2 on this level
	if isAMR2 then begin
	   if (amr1.levels[n].ngrids ne amr2.levels[n].ngrids) or $
	      (total(amr1.levels[n].idxlo ne amr2.levels[n].idxlo) ne 0) or $
	      (total(amr1.levels[n].idxhi ne amr2.levels[n].idxhi) ne 0) or $
	      (total(amr1.levels[n].gridspacing ne $
		   amr2.levels[n].gridspacing) ne 0) or $
	      (amr1.levels[n].nfab ne amr2.levels[n].nfab) then begin
		   print, 'Error: structures are incompatible on level ', $
			   strtrim(string(n),2)
		   amr_free, amr_new
		   return, -1
	   endif
	endif

	; create a new fab pointer for amr_new
	amr_new.levels[n].fabptr = ptr_new(*amr1.levels[n].fabptr)

	; for safety, null out all the data pointers within this fab pointer
	for m=0, amr_new.levels[n].nfab-1 do $
		(*amr_new.levels[n].fabptr)[m].dataptr = ptr_new()

	; now go through fabs
	for m=0, amr1.levels[n].nfab-1 do begin

		if isAMR2 then begin
		   ; do consistency checks on amr1 and amr2 in this fab
		   if (total( (*amr1.levels[n].fabptr)[m].idxlo ne $
			      (*amr2.levels[n].fabptr)[m].idxlo ) ne 0) or $
		      (total( (*amr1.levels[n].fabptr)[m].idxhi ne $
			      (*amr2.levels[n].fabptr)[m].idxhi ) ne 0) or $
		      (total( (*amr1.levels[n].fabptr)[m].idxtype ne $
			      (*amr2.levels[n].fabptr)[m].idxtype ) ne 0) or $
		      (total( (*amr1.levels[n].fabptr)[m].xlo ne $
			      (*amr2.levels[n].fabptr)[m].xlo ) ne 0) or $
		      (total( (*amr1.levels[n].fabptr)[m].xhi ne $
			      (*amr2.levels[n].fabptr)[m].xhi ) ne 0) $
		      then begin
			print, 'Error: structures are incompatible on', $
				' level ', strtrim(string(n),2), ', fab ', $
				strtrim(string(m),2)
			amr_free, amr_new
			return, -1
		   endif

		   ; form the sum, and set the new data pointer to point to it
		   datadiff = *(*amr1.levels[n].fabptr)[m].dataptr - $
			      *(*amr2.levels[n].fabptr)[m].dataptr
		endif else begin
		   datadiff = *(*amr1.levels[n].fabptr)[m].dataptr - amr2
		endelse

		; point to the result
		(*amr_new.levels[n].fabptr)[m].dataptr = ptr_new(datadiff)

	endfor
endfor

; return the new amr object
return, amr_new
end

