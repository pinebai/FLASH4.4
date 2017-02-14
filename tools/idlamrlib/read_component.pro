function read_component, componentName, level, plotDescriptor, $
	idxmin=idxmin, idxmax=idxmax, simplerefine = simplerefine, $
	verbose = verbose

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


; this routine reads the component componentName, refines
; where necessary to put it all on level level, and returns the result.
; The plot header must already have been read into plotDescriptor.
; Coarsening is done by bi-linear interpolation by default, or by
; simple sampling if simplerefine is set. By default the entire problem
; domain is read. If idxmin or idxmax are set (they should be ndim-element
; vectors), they specify the lower and upper corners of the domain to be
; read. idxmin and idxmax are specified in index space on the specified
; level of refinement.

; get component number and check that it exists
component = where(plotDescriptor.quantities eq componentName)
component = component[0]
if component eq -1 then begin
	print, 'Component '+componentName+' does not exist.'
	print, 'Valid components are:'
	print, plotDescriptor.quantities
	return, -1
endif
if level gt plotDescriptor.maxlevel then begin
	print, 'Maximum level available is ' + plotDescripto.maxlevel
	return, -1
endif

; figure out which family of multifab files contains this component,
; and which component number within that family it is
fileNameIndex = 0
while component ge $
	plotDescriptor.levels[level].nfilecomp[fileNameIndex] do begin
	component = component - $
		plotDescriptor.levels[level].nfilecomp[fileNameIndex]
	fileNameIndex = fileNameIndex + 1
endwhile

; set lower and upper index limits on quantities to be read
if not keyword_set(idxmin) then $
	idxmin = plotDescriptor.levels[level].idxlo
if not keyword_set(idxmax) then $
	idxmax = plotDescriptor.levels[level].idxhi
datasize = idxmax - idxmin + 1

; create array to hold result
if plotDescriptor.ndim eq 1 then $
	data=dblarr(datasize[0]) $
else if plotDescriptor.ndim eq 2 then $
	data=dblarr(datasize[0], datasize[1]) $
else if plotDescriptor.ndim eq 3 then $
	data=dblarr(datasize[0], datasize[1], datasize[2]) $
else begin
	print, 'Error: number of dimensions is not 1, 2, or 3.'
	return, -1
endelse

; now proceed through the levels, starting at the coarsest one and going
; until the target level is reached
for n=0, level do begin

	if keyword_set(verbose) then print, 'Reading level ', n, '...'

	; figure out the refinement ratio between this level and the
	; target level
	refineFactor = 1
	for m = n, level-1 do $
		refineFactor = refineFactor * plotDescriptor.refratio[m]

	; read header for the correct family of multifabs
	multifabHeader = read_multifabheader(n, fileNameIndex, $
		plotDescriptor)

	; go through all fabs associated with this header
	for m=0, multifabHeader.nfab-1 do begin

		; refine fab indices to target level
		fablo = multifabHeader.fabs[m].idxlo * refineFactor
		fabhi = (multifabHeader.fabs[m].idxhi+1) * $
			refineFactor - 1
		fabsize = fabhi - fablo + 1

		; determine if this fab is within the domain we're
		; trying to retrieve
		if total( ((fablo ge idxmin) and (fablo le idxmax)) or $
			  ((fabhi ge idxmin) and (fabhi le idxmax)) $
			) eq plotDescriptor.ndim $
		then begin

			; we do need to read this fab
			fab = read_fabcomponent(m, component, $
				plotDescriptor, multifabHeader)

			; refine fab to target level
			fabsize = fabhi - fablo + 1
			if keyword_set(simplerefine) then begin
				if plotDescriptor.ndim eq 1 then $
					fabref = rebin(fab, fabsize[0], $
						/sample) $
				else if plotDescriptor.ndim eq 2 then $
					fabref = rebin(fab, fabsize[0], $
						fabsize[1], /sample) $
				else if plotDescriptor.ndim eq 3 then $
					fabref = rebin(fab, fabsize[0], $
					fabsize[1], fabsize[2], /sample)
			endif else begin
				if plotDescriptor.ndim eq 1 then $
					fabref = rebin(fab, fabsize[0]) $
				else if plotDescriptor.ndim eq 2 then $
					fabref = rebin(fab, fabsize[0], $
					fabsize[1]) $
				else if plotDescriptor.ndim eq 3 then $
					fabref = rebin(fab, fabsize[0], $
					fabsize[1], fabsize[2])
			endelse

			; construct max and min within part of data
			; box that overlaps this fab
			datamin = fablo - idxmin
			datamin = (datamin ge 0) * datamin
			datamax = fabhi - idxmin
			datamax = (datamax lt datasize) * datamax + $
				(datamax ge datasize) * (datasize - 1)

			; construct max and min within this fab that
			; overlaps the data box
			fabmin = idxmin - fablo
			fabmin = (fabmin ge 0) * fabmin
			fabmax = idxmax - fablo
			fabmax = (fabmax lt fabsize) * fabmax + $
				(fabmax ge fabsize) * (fabsize - 1)

			; insert fab into data array
			if plotDescriptor.ndim eq 1 then $
				data[datamin[0] : datamax[0]] = $
				fabref[fabmin[0] : fabmax[0]] $
			else if plotDescriptor.ndim eq 2 then $
				data[datamin[0] : datamax[0], $
				     datamin[1] : datamax[1]] = $
				fabref[fabmin[0] : fabmax[0], $
				       fabmin[1] : fabmax[1]] $
			else if plotDescriptor.ndim eq 3 then $
				data[datamin[0] : datamax[0], $
				     datamin[1] : datamax[1], $
				     datamin[2] : datamax[2]] = $
				fabref[fabmin[0] : fabmax[0], $
				       fabmin[1] : fabmax[1], $
				       fabmin[2] : fabmax[2]]
		endif

	endfor

endfor

return, data

end