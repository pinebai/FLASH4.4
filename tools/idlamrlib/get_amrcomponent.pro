function get_amrcomponent, dirname, componentName

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

; routine to read a single amr component from a plot file / directory
; and store it in an IDL structure. This structure can then be passed 
; to other amrlib routines to do various useful things. The argument
; dirname is the name of the plot directory to be read, and componetName
; is the name of the component to be read.

; read plotfile header
plotdescriptor = read_amrheader(dirname)
if (plotdescriptor.ndim ne 1) and (plotdescriptor.ndim ne 2) and $
	(plotdescriptor.ndim ne 3) then begin
	print, 'Error: only 1, 2, or 3 dimensional plots are supported.'
	return, -1
endif
ndim = plotdescriptor.ndim

; get component number and check that it exists
component = where(plotDescriptor.quantities eq componentName)
component = component[0]
if component eq -1 then begin
	print, 'Component '+componentName+' does not exist.'
	print, 'Valid components are:'
	print, plotDescriptor.quantities
	return, -1
endif

; figure out which family of multifab files contains this component,
; and which component number within that family it is
fileNameIndex = 0
while component ge $
	plotDescriptor.levels[0].nfilecomp[fileNameIndex] do begin
	component = component - $
		plotDescriptor.levels[0].nfilecomp[fileNameIndex]
	fileNameIndex = fileNameIndex + 1
endwhile

; read number of sink particles if sink file exists
nsink=0
if strpos(dirname, '/', /reverse_search) eq strlen(dirname)-1 then $
	sinkname=dirname+'SinkParticles' $
else sinkname=dirname+'/SinkParticles'
if (file_test(sinkname)) then begin
	openr, fp, sinkname, /get_lun
	readf, fp, nsink
	free_lun, fp
endif

; read number of star particles if star file exists
nstar=0
if strpos(dirname, '/', /reverse_search) eq strlen(dirname)-1 then $
	starname=dirname+'StarParticles' $
else starname=dirname+'/StarParticles'
if (file_test(starname)) then begin
	openr, fp, starname, /get_lun
	readf, fp, nstar
	free_lun, fp
endif

; create a new structure in which to hold the amr data

; establish a structure to hold sink particles
sinkParticleDescriptor = { x:dblarr(ndim), p:dblarr(ndim), $
	j:dblarr(ndim), m:double(0.0) }

; establish a structure to hold star particles
starParticleDescriptor = { m:double(0.0), x:dblarr(ndim), p:dblarr(ndim), $
	j:dblarr(ndim), mlast:double(0.0), r:double(0.0), mdeut:double(0.0), $
	n:double(0.0), mdot:double(0.0), burnState:0 }

; the lowest level structure is a fab. It holds a single grid with a
; specified index range, and physical range, and data.

fabDescriptor = { idxlo:lonarr(ndim), idxhi:lonarr(ndim), $
	idxtype:lonarr(ndim), xlo:dblarr(ndim), xhi:dblarr(ndim), $
	dataptr:ptr_new(/allocate_heap) }
ptr_free, fabDescriptor.dataptr

; note: dataptr:ptr_new(/allocate_heap) is roughly equivalent to void *dataptr
; However, IDL insists on allocating valid memory to the pointer, which we must
; free up since we don't want it and will do our own memory assignment later.

; the next level object is a level. It consists of a series of fabs. It
; also stores various information about that level.

levelDescriptor = { level:0, ngrids:0, leveltime:0.0d, levelsteps:0L,  $
	idxlo:lonarr(ndim), idxhi:lonarr(ndim), $
	periodicity:lonarr(ndim), gridspacing:dblarr(ndim), $
	grids:plotdescriptor.levels[0].grids, $
	nfab:0, fabptr:ptr_new(/allocate_heap) }
ptr_free, levelDescriptor.fabptr	; free up dummy memory

; The highest level object is an amr descriptor. It consistents of a
; series of levels, and also stores global information about the amr grid.

; We have to have four cases here because IDL doesn't have a goddamn void *
; pointer!

if (nsink eq 0) then begin
	if (nstar eq 0) then begin
		; no stars or sinks
		amrDescriptor = { name:plotdescriptor.name, $
		componentName:componentName, $
		version:plotdescriptor.version, $
		ndim:plotdescriptor.ndim, $
		time:plotdescriptor.time, $
		maxlevel:plotdescriptor.maxlevel, $
		boxmin:plotdescriptor.boxmin, boxmax:plotdescriptor.boxmax, $
		refratio:plotdescriptor.refratio, $
		idxlo:plotdescriptor.idxlo, idxhi:plotdescriptor.idxhi, $
		idxtype:intarr(plotdescriptor.ndim), $
		periodicity:plotdescriptor.periodicity, $
		levelsteps:plotdescriptor.levelsteps, $
		gridspacing:plotdescriptor.gridspacing, $
		coordtype:plotdescriptor.coordtype, $
		levels:replicate(levelDescriptor,plotdescriptor.maxlevel+1), $
		nsink:0, sinkparticles:ptr_new(), $
		nstar:0, starparticles:ptr_new() }
	endif else begin
		; stars but no sinks
		amrDescriptor = { name:plotdescriptor.name, $
		componentName:componentName, $
		version:plotdescriptor.version, $
		ndim:plotdescriptor.ndim, $
		time:plotdescriptor.time, $
		maxlevel:plotdescriptor.maxlevel, $
		boxmin:plotdescriptor.boxmin, boxmax:plotdescriptor.boxmax, $
		refratio:plotdescriptor.refratio, $
		idxlo:plotdescriptor.idxlo, idxhi:plotdescriptor.idxhi, $
		idxtype:intarr(plotdescriptor.ndim), $
		periodicity:plotdescriptor.periodicity, $
		levelsteps:plotdescriptor.levelsteps, $
		gridspacing:plotdescriptor.gridspacing, $
		coordtype:plotdescriptor.coordtype, $
		levels:replicate(levelDescriptor,plotdescriptor.maxlevel+1), $
		nsink:0, sinkparticles:ptr_new(), $
		nstar:nstar, $
		starparticles:replicate(starParticleDescriptor,nstar) }
	endelse
endif else begin
	if (nstar eq 0) then begin
		; sinks but no stars
		amrDescriptor = { name:plotdescriptor.name, $
		componentName:componentName, $
		version:plotdescriptor.version, $
		ndim:plotdescriptor.ndim, $
		time:plotdescriptor.time, $
		maxlevel:plotdescriptor.maxlevel, $
		boxmin:plotdescriptor.boxmin, boxmax:plotdescriptor.boxmax, $
		refratio:plotdescriptor.refratio, $
		idxlo:plotdescriptor.idxlo, idxhi:plotdescriptor.idxhi, $
		idxtype:intarr(plotdescriptor.ndim), $
		periodicity:plotdescriptor.periodicity, $
		levelsteps:plotdescriptor.levelsteps, $
		gridspacing:plotdescriptor.gridspacing, $
		coordtype:plotdescriptor.coordtype, $
		levels:replicate(levelDescriptor,plotdescriptor.maxlevel+1), $
		nsink:nsink, $
		sinkparticles:replicate(sinkParticleDescriptor,nsink), $
		nstar:0, starparticles:ptr_new() }
	endif else begin
		; both sinks and stars -- this probably shouldn't ever
		; happen, but we include the possibility anyway
		amrDescriptor = { name:plotdescriptor.name, $
		componentName:componentName, $
		version:plotdescriptor.version, $
		ndim:plotdescriptor.ndim, $
		time:plotdescriptor.time, $
		maxlevel:plotdescriptor.maxlevel, $
		boxmin:plotdescriptor.boxmin, boxmax:plotdescriptor.boxmax, $
		refratio:plotdescriptor.refratio, $
		idxlo:plotdescriptor.idxlo, idxhi:plotdescriptor.idxhi, $
		idxtype:intarr(plotdescriptor.ndim), $
		periodicity:plotdescriptor.periodicity, $
		levelsteps:plotdescriptor.levelsteps, $
		gridspacing:plotdescriptor.gridspacing, $
		coordtype:plotdescriptor.coordtype, $
		levels:replicate(levelDescriptor,plotdescriptor.maxlevel+1), $
		nsink:nsink, $
		sinkparticles:replicate(sinkParticleDescriptor,nstar), $
		nstar:nstar, $
		starparticles:replicate(starParticleDescriptor,nstar) }
	endelse
endelse

; now read all the data

for n=0, amrDescriptor.maxlevel do begin

	; read header file
	mfhdr = read_multifabHeader(n, fileNameIndex, plotdescriptor)

	; store useful info
	amrDescriptor.levels[n].level=n
	amrDescriptor.levels[n].ngrids = $
		plotdescriptor.levels[n].ngrids
	amrDescriptor.levels[n].leveltime = $
		plotdescriptor.levels[n].leveltime
	amrDescriptor.levels[n].levelsteps = $
		plotdescriptor.levels[n].levelsteps
	amrDescriptor.levels[n].idxlo = $
		plotdescriptor.levels[n].idxlo
	amrDescriptor.levels[n].idxhi = $
		plotdescriptor.levels[n].idxhi
	amrDescriptor.levels[n].periodicity = $
		plotdescriptor.levels[n].periodicity
	amrDescriptor.levels[n].gridspacing = $
		plotdescriptor.levels[n].gridspacing
	amrDescriptor.levels[n].grids = $
		plotdescriptor.levels[n].grids
	amrDescriptor.levels[n].nfab = mfhdr.nfab

	; we now know the number of fabs, so create an array of
	; fab descriptors for the fabptr to point to.
	amrDescriptor.levels[n].fabptr = $
		ptr_new(replicate(fabDescriptor,mfhdr.nfab))

	; loop through the fabs
	for l=0, amrDescriptor.levels[n].nfab-1 do begin

		; record properties of each fab. idxlo and idxhi
		; are the index limits, idxtype is the index type
		; (0 = cell centered, 1 = edge centered), and
		; xlo and xhi are the physical limits
		(*amrDescriptor.levels[n].fabptr)[l].idxlo = $
			mfhdr.fabs[l].idxlo
		(*amrDescriptor.levels[n].fabptr)[l].idxhi = $
			mfhdr.fabs[l].idxhi
		(*amrDescriptor.levels[n].fabptr)[l].idxtype = $
			mfhdr.fabs[l].idxtype
		(*amrDescriptor.levels[n].fabptr)[l].xlo = $
			plotdescriptor.boxmin + $
			plotdescriptor.gridspacing[*,n] * mfhdr.fabs[l].idxlo
		(*amrDescriptor.levels[n].fabptr)[l].xhi = $
			plotdescriptor.boxmin + $
			plotdescriptor.gridspacing[*,n] * $
			(mfhdr.fabs[l].idxhi + $
				(mfhdr.fabs[l].idxtype eq 0))

		; read in fab data, simultaneously allocating memory to hold it
		(*amrDescriptor.levels[n].fabptr)[l].dataptr = ptr_new( $
			read_fabcomponent(l, component, plotdescriptor, mfhdr))

	endfor

endfor

; Set the index type at the top level. In the plot file format, this
; is stored fab by fab, since different components can have different
; types (cell vs. edge centered) within the same plot file. However,
; since amrDescriptor applies to only a single component, it makes more
; sense for the index type to be stored at the top level -- it will be
; the same for every fab in any case.
amrDescriptor.idxtype = (*amrDescriptor.levels[0].fabptr)[0].idxtype

; read the sink particle data
if (nsink ne 0) then begin
	openr, fp, sinkname, /get_lun
	readf, fp, nsink
	sinkdata=dblarr(3*ndim+1)
	sinkptr=0
	for sinkptr=0, nsink-1 do begin
		readf, fp, sinkdata
		amrDescriptor.sinkparticles[sinkptr].m = sinkdata[0]
		amrDescriptor.sinkparticles[sinkptr].x = sinkdata[1:ndim]
		amrDescriptor.sinkparticles[sinkptr].p = $
			sinkdata[ndim+1:2*ndim]
		amrDescriptor.sinkparticles[sinkptr].j = $
			sinkdata[2*ndim+1:3*ndim]
	endfor
	free_lun, fp
endif

; read the star particle data
if (nstar ne 0) then begin
	openr, fp, starname, /get_lun
	readf, fp, nstar
	stardata=dblarr(3*ndim+7)
	starptr=0
	for sinkptr=0, nstar-1 do begin
		readf, fp, stardata
		amrDescriptor.starparticles[sinkptr].m = stardata[0]
		amrDescriptor.starparticles[sinkptr].x = stardata[1:ndim]
		amrDescriptor.starparticles[sinkptr].p = $
			stardata[ndim+1:2*ndim]
		amrDescriptor.starparticles[sinkptr].j = $
			stardata[2*ndim+1:3*ndim]
		amrDescriptor.starparticles[sinkptr].mlast = stardata[3*ndim+1]
		amrDescriptor.starparticles[sinkptr].r = stardata[3*ndim+2]
		amrDescriptor.starparticles[sinkptr].mdeut = stardata[3*ndim+3]
		amrDescriptor.starparticles[sinkptr].n = stardata[3*ndim+4]
		amrDescriptor.starparticles[sinkptr].mdot = stardata[3*ndim+5]
		amrDescriptor.starparticles[sinkptr].burnstate = $
			stardata[3*ndim+6]
	endfor
	free_lun, fp
endif

; return data
return, amrDescriptor
end
