function read_amrheader, dirname, verbose=verbose

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


; function to read the header file for and AMR plot directory

; open header file
if strpos(dirname, '/', /reverse_search) eq strlen(dirname)-1 then $
	hdrname=dirname+'Header' else hdrname=dirname+'/Header'
openr, fp, hdrname, /get_lun

; read header info about entire problem

version=''
readf, fp, version		; version number

nplot=0				; number of plotted quantites
readf, fp, nplot

quantities = strarr(nplot)
temp=''
for n=0, nplot-1 do begin	
	readf, fp, temp		; plotted quantity names
	quantities[n] = temp
endfor

ndim=0
readf, fp, ndim			; number of spatial dimensions

time=0.0
readf, fp, time			; time of plot file

maxlevel=0
readf, fp, maxlevel		; maximum level of refinement

boxmin=fltarr(ndim)
boxmax=fltarr(ndim)
readf, fp, boxmin		; positions of box corners
readf, fp, boxmax

readf, fp, temp			; refinement ratios
if maxlevel gt 0 then refratio=lonarr(maxlevel) else refratio=0
for n=0, maxlevel-1 do begin
	temp=strmid(temp, strpos(temp, "(")+1)
	temp1=strmid(temp, 0, strpos(temp, ","))
	refratio[n]=fix(temp1)
endfor

readf, fp, temp			; problem domain in index space
idxlo=lonarr(ndim,maxlevel+1)
idxhi=lonarr(ndim,maxlevel+1)
periodicity=intarr(ndim,maxlevel+1)
for n=0, maxlevel do begin
	temp=strmid(temp, strpos(temp, "((")+2)
	for l=0, ndim-2 do begin
		temp1=strmid(temp, 0, strpos(temp, ","))
		idxlo[l,n]=long(temp1)
		temp=strmid(temp, strpos(temp, ",")+1)
	endfor
	temp1=strmid(temp, 0, strpos(temp, ")"))
	idxlo[l,n]=long(temp1)	
	temp=strmid(temp, strpos(temp, "(")+1)	
	for l=0, ndim-2 do begin
		temp1=strmid(temp, 0, strpos(temp, ","))
		idxhi[l,n]=long(temp1)
		temp=strmid(temp, strpos(temp, ",")+1)
	endfor
	temp1=strmid(temp, 0, strpos(temp, ")"))
	idxhi[l,n]=long(temp1)	
	temp=strmid(temp, strpos(temp, "(")+1)	
	for l=0, ndim-2 do begin
		temp1=strmid(temp, 0, strpos(temp, ","))
		periodicity[l,n]=long(temp1)
		temp=strmid(temp, strpos(temp, ",")+1)
	endfor
	temp1=strmid(temp, 0, strpos(temp, ")"))
	periodicity[l,n]=long(temp1)	
endfor

levelsteps=lonarr(maxlevel+1)	; number of steps on each level
readf, fp, levelsteps

dx=fltarr(ndim,maxlevel+1)	; grid spacing on each level
readf, fp, dx

coordtype=0
readf, fp, coordtype		; coordinate system type; 0 = cartesian

readf, fp, temp			; dummy line


; First step done we've read all the info about the overall structure as
; opposed to invidual levels

; if verbose dump output
if keyword_set(verbose) then begin
	print, 'File format version: ', version
	print, 'Plotted quantities:'
	for n=0, nplot-1 do print, '  ', quantities[n]
	print, 'Number of dimensions: ', strtrim(string(ndim), 2)
	print, 'Time: ', strtrim(string(time),2)
	print, 'Maximum level of refinement: ', strtrim(string(maxlevel),2)
	print, 'Box minima: ', boxmin
	print, 'Box maxima: ', boxmax
	if maxlevel gt 0 then print, 'Refinement ratios: ', refratio
	for n=0,maxlevel do begin
		print, 'Level ', strtrim(string(n),2),':'
		print, '  Domain min = ', reform(idxlo[*,n])
		print, '  Domain max = ', reform(idxhi[*,n])
		print, '  Domain periodicity = ', reform(periodicity[*,n])
		print, '  Number of steps = ', $
			strtrim(string(levelsteps[n]),2)
	endfor
	print, 'Coordinate system type: ', strtrim(string(coordtype),2)
endif

; Now prepare to read about individual levels and grids

; create a generic grid descriptor
gridDescriptor = { xlo:dblarr(ndim), xhi:dblarr(ndim) }

; Loop over levels. The first time we do this, we just have to figure
; out the maximum number of grids on any level. This is because
; idl won't allow structures to contain arrays of unspecified
; size, so we need to get the largest number of grids that will be needed
; on any level, THEN define the level descriptor, THEN define the
; plot file descriptor. This could be accomplished with pointers in
; IDL, but IDL pointers are clumsy and do not allow the same nifty
; array operations that IDL arrays allow. Therefore we use this
; bit of messy coding here to avoid more messy code later on.

; save file position
posn=0L
point_lun, -fp, posn

; get maximum number of grids
maxgrids=0
for n=0, maxlevel do begin

	level=-1
	grids=0
	levtime=0.0
	readf, fp, level, grids, levtime	; read current level

	if (grids gt maxgrids) then maxgrids = grids

	readf, fp, temp		; number of steps on this level, already stored

	gridboundaries = dblarr(2,ndim)
	for m=0, grids-1 do begin		; read data on each grid
		readf, fp, gridboundaries
	endfor

	ncomp=0
	nfiles=0
	repeat begin				; read filenames
		readf, fp, temp
		nfiles = nfiles + 1

		; we need to open this multifab header file and read
		; the first 3 lines to find out how many components
		; are associated with it
		if strpos(dirname, '/', /reverse_search) $
			eq strlen(dirname)-1 then $
			multifabname=dirname+temp+'_H' else $
			multifabname=dirname+'/'+temp+'_H'
		openr, fp1, multifabname, /get_lun
		dummy = 0
		nmultifabcomp = 0
		readf, fp1, dummy
		readf, fp1, dummy
		readf, fp1, nmultifabcomp
		free_lun, fp1

		; figure out how many total components there are in
		; all the multifab header files we've read already,
		; and keep reading only if there are quantities we
		; haven't read yet
		ncomp = ncomp + nmultifabcomp
		if ncomp gt nplot then begin
			print, 'Error: inconsistent number of ' + $
				'quantities in multifab header files.'
			return, -1
		endif
	endrep until ncomp eq nplot

endfor

; ok, now we can create a level descriptor
levelDescriptor = { level:0, ngrids:0, leveltime:0.0d, levelsteps:0L,  $
	idxlo:lonarr(ndim), idxhi:lonarr(ndim), $
	periodicity:lonarr(ndim), gridspacing:dblarr(ndim), $
	grids:replicate(gridDescriptor,maxgrids), $
	nfiles:nfiles, filenames:strarr(nfiles), $
	nfilecomp:intarr(nfiles) } 

; create the plot descriptor
thisPlot = { name:dirname, version:version, nplot:nplot, $
	quantities:quantities, ndim:ndim, time:time, $
	maxlevel:maxlevel, $
	boxmin:boxmin, boxmax:boxmax, refratio:refratio, $
	idxlo:idxlo, idxhi:idxhi, $
	periodicity:periodicity, levelsteps:levelsteps, $
	gridspacing:dx, coordtype:coordtype, $
	levels:replicate(levelDescriptor,maxlevel+1) }


; now loop through levels a second time and actually store the data

; go to previous file position
point_lun, fp, posn

for n=0, maxlevel do begin

	readf, fp, level, grids, levtime	; read level, grids, time
	thisPlot.levels[level].level = level	; store
	thisPlot.levels[level].ngrids = grids
	thisPlot.levels[level].leveltime = levtime

	levelsteps=0L
	readf, fp, levelsteps	; number of steps on this level
	thisPlot.levels[level].levelsteps = levelsteps

	; load the stuff already recorded into this level descriptor
	thisPlot.levels[level].idxlo = $
		reform(thisPlot.idxlo[*,level])
	thisPlot.levels[level].idxhi = $
		reform(thisPlot.idxhi[*,level])
	thisPlot.levels[level].periodicity = $
		reform(thisPlot.periodicity[*,level])
	thisPlot.levels[level].gridspacing = $
		reform(thisPlot.gridspacing[*,level])

	; now read data about individual grids
	for m=0, grids-1 do begin		; read data on each grid
		readf, fp, gridboundaries
		thisPlot.levels[level].grids[m].xlo = $
			reform(gridboundaries[0,*])
		thisPlot.levels[level].grids[m].xhi = $
			reform(gridboundaries[1,*])
	endfor

	; read file names for this level
	fname=''
	for m=0, nfiles-1 do begin
		readf, fp, fname	; read Level_N/...
		thisPlot.levels[level].filenames[m] = fname

		; again, open the file to read the number of components
		if strpos(dirname, '/', /reverse_search) $
			eq strlen(dirname)-1 then $
			multifabname=dirname+fname+'_H' else $
			multifabname=dirname+'/'+fname+'_H'
		openr, fp1, multifabname, /get_lun
		readf, fp1, dummy
		readf, fp1, dummy
		readf, fp1, nmultifabcomp
		free_lun, fp1

		; store number of components of this file
		thisPlot.levels[level].nfilecomp[m] = nmultifabcomp

	endfor

endfor


; close header file
free_lun, fp

return, thisPlot

end
