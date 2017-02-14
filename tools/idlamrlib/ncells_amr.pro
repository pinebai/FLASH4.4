function ncells_amr, amr, maxlevel=maxlevel

; Copyright Mark Krumholz (2001)
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


; this routine finds the number of cells on each level in an AMR grid

if not keyword_set(maxlevel) then maxlevel=amr.maxlevel
ncells=lonarr(maxlevel+1)

; go through amr structure
for n=0, maxlevel do begin
	for m=0, amr.levels[n].nfab-1 do begin

		; add number of cells to total
		ncells[n] = ncells[n] + $
			n_elements(*(*amr.levels[n].fabptr)[m].dataptr)

	endfor
endfor

return, ncells

end



