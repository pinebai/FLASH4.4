pro raster, z, x, y, anisotropic=anisotropic, $
	xrange=xrange, yrange=yrange, nocolorbar=nocolorbar, $
	verbose=verbose, zrange=zrange, localrange=localrange, $
	snaptogrid=snaptogrid, _extra = extra

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


; this routine produces a raster image in the plot window.
; Limitations: the routine assumes that z is uniformly spaced and x
; 	and y are uniformly spaced and strictly increasing. It uses
;	x and y only to label the axes, not to determine where pixels
;	are placed in the raster image.

; read arguments
zref=reform(z)
sz=size(zref)
if sz[0] ne 2 then begin
	print, 'Error: image must be a two dimensional array'
	return
endif
if n_params() eq 1 then begin
	x = lindgen(sz[1])
	y = lindgen(sz[2])
endif else if n_params() ne 3 then begin
	print, 'Error: raster requires 1 or 3 parameters'
	return
endif

; set up x and y min, max, and spacing
ncells=sz[1:2]
gridspacing = [(x[ncells[0]-1]-x[0])/(ncells[0]-1), $
	       (y[ncells[1]-1]-y[0])/(ncells[1]-1)]
xmin=x[0] - gridspacing[0]/2.0
xmax=x[ncells[0]-1] + gridspacing[0]/2.0
ymin=y[0] - gridspacing[1]/2.0
ymax=y[ncells[1]-1] + gridspacing[1]/2.0
boxmin=[xmin, ymin]
boxmax=[xmax, ymax]

; get keywords
if not keyword_set(noerase) then noerase=0
if not keyword_set(anisotropic) then isotropic=1 else isotropic=0
if not keyword_set(xrange) then xlim=[xmin, xmax] else xlim=xrange
if not keyword_set(yrange) then ylim=[ymin, ymax] else ylim=yrange

; if asked to do so, snap range to grid
if keyword_set(xrange) and keyword_set(snaptogrid) then begin
	xlim[0] = floor( (xlim[0] - xmin) / gridspacing[0] ) * $
		  gridspacing[0] + xmin
	xlim[1] = ceil( (xlim[1] - xmin) / gridspacing[0] ) * $
		  gridspacing[0] + xmin
	if keyword_set(verbose) then print, 'New x range = ', xlim
endif
if keyword_set(yrange) and keyword_set(snaptogrid) then begin
	ylim[0] = floor( (ylim[0] - ymin) / gridspacing[1] ) * $
		  gridspacing[1] + ymin
	ylim[1] = ceil( (ylim[1] - ymin) / gridspacing[1] ) * $
		  gridspacing[1] + ymin
	if keyword_set(verbose) then print, 'New y range: ', ylim
endif

; get some information about the graphics environment
bangmulti = !p.multi
tempmulti = !p.multi
tempmulti[1:3] = tempmulti[1:3] + (tempmulti[1:3] eq 0)
xmargindev = !x.margin * !d.x_ch_size
ymargindev = !y.margin * !d.y_ch_size

; Convert data to normal coordinates. If we're on a device that uses
; windows, make sure a valid window is open before attempting this.
devname = !d.name
devname = strmid(devname, 0, 3)
devname = strupcase(devname)
case devname of
	'MAC': if !d.window eq -1 then window, 0
	'WIN': if !d.window eq -1 then window, 0
	'X': if !d.window eq -1 then window, 0
	ELSE:
endcase
marginnorm = convert_coord(xmargindev, ymargindev, /device, /to_normal)

; see how many colors are available
ncolors = !d.table_size

; draw axes
if keyword_set(nocolorbar) then $
	plot, [0,0], [0,0], /nodata, xrange=xlim, yrange=ylim, $
		isotropic=isotropic, _extra = extra $
else begin
	; in the case with a colorbar, we must deal with multiple plots
	; (!p.multi settings) manually

	; set the row and column for this plot
	if tempmulti[4] eq 0 then begin
		plotrow = tempmulti[2] - tempmulti[0] / tempmulti[1] - 1
		plotcol = tempmulti[0] mod tempmulti[1]
	endif else begin
		plotrow = tempmulti[2] - (tempmulti[0] mod tempmulti[2]) - 1
		plotcol = tempmulti[0] / tempmulti[2]
	endelse

	; set the plot box position
	plotpos = fltarr(4)
	plotpos[0] = marginnorm[0] + float(plotcol) / tempmulti[1]
	plotpos[1] = marginnorm[1] + float(plotrow) / tempmulti[2]
	plotpos[2] = (plotcol+0.7) / tempmulti[1] - marginnorm[2]
	plotpos[3] = float(plotrow+1) / tempmulti[2] - marginnorm[3]

	; draw the plot box
	plot, [0,0], [0,0], /nodata, xrange=xlim, yrange=ylim, $
		isotropic=isotropic, _extra = extra, $
		position = plotpos
endelse

; What we do here depends on whether we have a device with scalable
; pixels or not. If we have scalable pixels, we can just construct
; an image box where each cell of the amr grid on the finest level
; represents a single pixel. We then rely on the scaling of pixels
; to make that image correctly overlay the plot window. If we don't
; have scalable pixels, the number of cells per pixel will in be a
; complicated function of the size of the graphics device and the
; resolution of the amr grid. For reference, X windows does not have
; scalable pixels, postscript does.

scalablePixels = !d.flags mod 2

if scalablePixels then begin

	; we have scalable pixels

	; Figure out the max and min cell index corresponding to the
	; xlim and ylim limits we've been given.
	xidxlo = floor( (xlim[0] - xmin)/gridspacing[0] )
	xidxhi = floor( (xlim[1] - xmin)/gridspacing[0] ) - 1
	yidxlo = floor( (ylim[0] - ymin)/gridspacing[1] )
	yidxhi = floor( (ylim[1] - ymin)/gridspacing[1] ) - 1

	; Set high and low indices for the image box
	imglo = long([xidxlo, yidxlo])
	imghi = long([xidxhi, yidxhi])

endif else begin

	; we do not have scalable pixels

	; get max and min for image box in pixels
	imglo = convert_coord(xlim[0], ylim[0], /data, /to_device)
	imghi = convert_coord(xlim[1], ylim[1], /data, /to_device)
	imglo=round(imglo[0:1])
	imghi=round(imghi[0:1])
	if (imghi[0] - imglo[0]) lt 2 then begin
		print, 'Error: requested x axis range is too small.'
		return
	endif
	if (imghi[1] - imglo[1]) lt 2 then begin
		print, 'Error: requested y axis range is too small.'
		return
	endif

endelse

; construct the image box
img=fltarr(imghi[0]-imglo[0]+1, imghi[1]-imglo[1]+1)

; figure out which region of the data box overlaps with the image box

; first figure this out in data coordinates
datalo=fltarr(2)
datahi=fltarr(2)
datalo[0] = (xmin gt xlim[0]) * xmin + (xmin le xlim[0]) * xlim[0]
datalo[1] = (ymin gt ylim[0]) * ymin + (ymin le ylim[0]) * ylim[0]
datahi[0] = (xmax le xlim[1]) * xmax + (xmax gt xlim[1]) * xlim[1]
datahi[1] = (ymax le ylim[1]) * ymax + (ymax gt ylim[1]) * ylim[1]

; if there is no overlap between the image box and the data,
; then there is nothing to display so we return
if (datalo[0] gt xmax) or (datahi[0] lt xmin) or $
   (datalo[1] gt ymax) or (datahi[1] lt ymin) $
	then return

; convert limits in data coordinates to indices with the fab
idxlo = round( (datalo-boxmin)/gridspacing )
idxhi = round( (datahi-boxmin)/gridspacing ) - 1
idxlo = (idxlo ge 0) * idxlo
idxhi = (idxhi gt (ncells-1)) * (ncells-1) + $
	(idxhi le (ncells-1)) * idxhi

; We must convert limits in data coordinates to indices in
; the image array. How we do this depends on whether we
; have scalable pixels or not.

if scalablePixels then begin

	; we have scalable pixels
	imgdatalo = (datalo-boxmin)/(boxmax-boxmin) * ncells
	imgdatahi = (datahi-boxmin)/(boxmax-boxmin) * ncells - 1
	imgidxlo = long(imgdatalo-imglo)
	imgidxhi = long(imgdatahi-imglo)		

endif else begin

	; we do not have scalable pixels
	imgdatalo = convert_coord(datalo[0], datalo[1], $
		/data, /to_device)
	imgdatahi = convert_coord(datahi[0], datahi[1], $
		/data, /to_device)
	imgidxlo = round(imgdatalo-imglo)
	imgidxhi = round(imgdatahi-imglo)
	imgidxlo = (imgidxlo ne -1) * imgidxlo
	imgidxhi = (imgidxhi ne (imghi-imglo+1)) * imgidxhi + $
		   (imgidxhi eq (imghi-imglo+1)) * $
			(imghi-imglo)

endelse

; copy data into image box
img[imgidxlo[0]:imgidxhi[0], imgidxlo[1]:imgidxhi[1]] = $
	congrid( zref[idxlo[0]:idxhi[0], idxlo[1]:idxhi[1]], $
		 imgidxhi[0]-imgidxlo[0]+1, imgidxhi[1]-imgidxlo[1]+1)

; set the range to be used to assign colors
if n_elements(zrange) eq 0 then begin
	if not keyword_set(localrange) then begin
		maxrange = max(z)
		minrange = min(z)
		if keyword_set(verbose) then begin
			print, 'Global minimum = ', minrange
			print, 'Global maximum = ', maxrange
		endif
	endif else begin
		maxrange = max(img)
		minrange = min(img)
		if keyword_set(verbose) then begin
			print, 'Local minimum = ', minrange
			print, 'Local maximum = ', maxrange
		endif
	endelse
endif else begin
	maxrange = zrange[1]
	minrange = zrange[0]
	if keyword_set(verbose) then begin
		print, 'User-provided minimum = ', minrange
		print, 'User-provided maximum = ', maxrange
	endif
endelse

; draw the color bar
if not keyword_set(nocolorbar) then begin

	; set up color bar string formatting
	if (abs(maxrange) lt 1.0e4) and (abs(maxrange) gt 1.0e-4) then $
		format = '(f8.2)' $
	else $
		format = '(e10.2)'

	; set the colorbar position
	colorbarpos = fltarr(4)
	colorbarpos[0] = (plotcol + 0.88) / tempmulti[1]
	colorbarpos[1] = plotpos[1]
	colorbarpos[2] = (plotcol+0.95) / tempmulti[1]
	colorbarpos[3] = plotpos[3]

	; draw the color bar
	colorbar, /vertical, maxrange = maxrange, minrange = minrange, $
		format = format, ncolors=ncolors, position=colorbarpos
endif

; put pixels into specified range
if n_elements(zrange) ne 0 then begin
	img = (img lt minrange) * minrange + (img ge minrange) * img
	img = (img gt maxrange) * maxrange + (img le maxrange) * img
endif

; scale image to color table
; code to deal with number of image colors lifted from Fanning colorbar.pro
;if scalablePixels then begin
;	oldDevice = !D.NAME

	; What kind of computer are we using? SET_PLOT to appropriate
	; display device.
;	thisOS = !VERSION.OS_FAMILY
;	thisOS = STRMID(thisOS, 0, 3)
;	thisOS = STRUPCASE(thisOS)
;	CASE thisOS of
;		'MAC': SET_PLOT, thisOS
;		'WIN': SET_PLOT, thisOS
;		ELSE: SET_PLOT, 'X'
;	ENDCASE

	; Here is how many colors we should use.
;	ncolors = !D.TABLE_SIZE
;	SET_PLOT, oldDevice

;ENDIF ELSE ncolors = !D.TABLE_SIZE

; scale the image to the color table
img1 = fix((img - minrange) / (maxrange - minrange) * (ncolors-1))

; Turn off device decomposition for scalable pixels
if not scalablePixels then begin
	if !Version.Release ge 5.2 then device, get_decomposed=thisdecomposed
	device, decomposed=0
endif

; Now display the image. How we do this again depends on whether we
; have scalable pixels or not.

if scalablePixels then begin

	; we have scalable pixels

	; get box limits in normal coordinates
	imglonorm = convert_coord(xlim[0], ylim[0], /data, /to_norm)
	imghinorm = convert_coord(xlim[1], ylim[1], /data, /to_norm)

	; show image
	tv, img1, imglonorm[0], imglonorm[1], $
		xsize = imghinorm[0] - imglonorm[0], $
		ysize = imghinorm[1] - imglonorm[1], $
		/norm

endif else begin

	; we do not have scalable pixels
	tv, img1, imglo[0], imglo[1]
endelse

; turn decomposition back to previous setting
if not scalablePixels then begin
	if !Version.Release ge 5.2 then device, decomposed=thisdecomposed
endif

; overlay the plot axes again
!p.multi=bangmulti
if keyword_set(nocolorbar) then $
	plot, [0,0], [0,0], /noerase, /nodata, xrange=xlim, yrange=ylim, $
		isotropic=isotropic, _extra = extra $
else $
	plot, [0,0], [0,0], /noerase, /nodata, xrange=xlim, yrange=ylim, $
		isotropic=isotropic, _extra = extra, $
		position = plotpos

; increment !p.multi, since the plot call with /noerase set doesn't do
; that automatically
!p.multi[0] = !p.multi[0] + 1
if !p.multi[0] eq tempmulti[1]*tempmulti[2]*tempmulti[3] then !p.multi[0]=0

return
end

