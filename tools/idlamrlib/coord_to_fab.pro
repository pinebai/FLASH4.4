function coord_to_fab, coord, amr

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


; This function finds all fabs, and the indices within them, corresponding
; to a given physical position. The input is coord, an ndim element real
; vector, and the amr object to be search. The output is an array of
; [nlevel x (2*ndim + 2)] elements, where nlevel is the number of levels
; that have fabs corresponding to that position. The output is organized
; as follows:
;	level of fab	fab number	[x,y,z,...] indices within fab
;					[x,y,z,...] indices within level
; So, for example, a return value of:
;	0	4	15	12	11
;	1	31	11	8	6
; means the point corresponds to cell [15,12,11] of the 4th fab on
; level 0 and cell [11,8,6] of the 31st fab on level 1. If the specified
; point is outside the amr domain, the routine returns -1.

; initialize results pointer
reslevels = 0

; go through amr structure
for n=0, amr.maxlevel do begin
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
			offset = coord - (*amr.levels[n].fabptr)[m].xlo
			fabidx = round((offset - 0.5*amr.gridspacing[*,n]) / $
				    amr.gridspacing[*,n])
			fabidx = fabidx * (fabidx ne -1)
			fabidxmax = (*amr.levels[n].fabptr)[m].idxhi - $
				(*amr.levels[n].fabptr)[m].idxlo
			fabidx = fabidx * (fabidx ne fabidxmax) + $
				fabidxmax * (fabidx eq fabidxmax)
			reslevels = reslevels + 1
			resnew = intarr(2*amr.ndim+2, reslevels)	
			if reslevels ne 1 then resnew[*, 0:reslevels-2] = res
			resnew[0, reslevels-1] = n
			resnew[1, reslevels-1] = m
			if amr.ndim eq 1 then begin
				resnew[2, reslevels-1] = fabidx[0]
				resnew[3, reslevels-1] = fabidx[0] + $
					(*amr.levels[n].fabptr)[m].idxlo[0]
			endif else if amr.ndim eq 2 then begin
				resnew[2, reslevels-1] = fabidx[0]
				resnew[3, reslevels-1] = fabidx[1]
				resnew[4, reslevels-1] = fabidx[0] + $
					(*amr.levels[n].fabptr)[m].idxlo[0]
				resnew[5, reslevels-1] = fabidx[1] + $
					(*amr.levels[n].fabptr)[m].idxlo[1]
			endif else begin
				resnew[2, reslevels-1] = fabidx[0]
				resnew[3, reslevels-1] = fabidx[1]
				resnew[4, reslevels-1] = fabidx[2]
				resnew[5, reslevels-1] = fabidx[0] + $
					(*amr.levels[n].fabptr)[m].idxlo[0]
				resnew[6, reslevels-1] = fabidx[1] + $
					(*amr.levels[n].fabptr)[m].idxlo[1]
				resnew[7, reslevels-1] = fabidx[2] + $
					(*amr.levels[n].fabptr)[m].idxlo[2]
			endelse
			res = resnew
			m = amr.levels[n].nfab
		endelse
	endfor
endfor

if reslevels ne 0 then return, res else return, -1

end
		
				



