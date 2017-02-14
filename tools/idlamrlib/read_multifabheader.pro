function read_multifabheader, level, filenameIdx, plotDescriptor

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


; function to read the header (_H) file associated with a set of
; multifab files on level level corresponding to filename number
; filenameIdx.

ndim=plotDescriptor.ndim

; open file
basename=plotDescriptor.levels[level].filenames[filenameIdx]
fname=plotDescriptor.name + "/" + basename + "_H"
openr, fp, fname, /get_lun

version=0		; output version number
readf, fp, version

method=0		; output method, with 0 = one file per CPU
readf, fp, method	; right now only method 0 works

ncomp=0			; number of components in these multifab files
readf, fp, ncomp

ngrow=0			; number of "ghost cells" included
readf, fp, ngrow	; right now only 0 works

temp=''
readf, fp, temp		; lines of text that must be parsed
temp=strmid(temp, strpos(temp, "(")+1)	; remove leading parenthesis
temp1=strmid(temp, 0, strpos(temp, " "))	; select portion before space
nfab=long(temp1)		; number of boxes

temp=strmid(temp, strpos(temp, " ")+1)	; select portion after space
hash_sig=long(temp)	; I have no idea what this is

; now we know enough to build the descriptors

; create a generic descriptor for a box
fabDescriptor = { idxlo:lonarr(plotdescriptor.ndim), $
	idxhi:lonarr(plotdescriptor.ndim), $
	idxtype:intarr(plotdescriptor.ndim), $
	filename:'', offset:0L, min:dblarr(ncomp), max:dblarr(ncomp) }

; create a descriptor for this multifab header
mfheaderDescriptor = { fullname:fname, basename:basename, $
	version:version, method:method, $
	ncomp:ncomp, ngrow:ngrow, nfab:nfab, hash_sig:hash_sig, $
	fabs:replicate(fabDescriptor, nfab) }

; read through the boxes
for n=0, nfab-1 do begin
	readf, fp, temp			; read line of text
	temp=strmid(temp, strpos(temp, "((")+2)
	for l=0, ndim-2 do begin
		temp1=strmid(temp, 0, strpos(temp, ","))
		mfheaderDescriptor.fabs[n].idxlo[l]=long(temp1)
		temp=strmid(temp, strpos(temp, ",")+1)
	endfor
	temp1=strmid(temp, 0, strpos(temp, ")"))
	mfheaderDescriptor.fabs[n].idxlo[l]=long(temp1)	
	temp=strmid(temp, strpos(temp, "(")+1)	
	for l=0, ndim-2 do begin
		temp1=strmid(temp, 0, strpos(temp, ","))
		mfheaderDescriptor.fabs[n].idxhi[l]=long(temp1)
		temp=strmid(temp, strpos(temp, ",")+1)
	endfor
	temp1=strmid(temp, 0, strpos(temp, ")"))
	mfheaderDescriptor.fabs[n].idxhi[l]=long(temp1)	
	temp=strmid(temp, strpos(temp, "(")+1)	
	for l=0, ndim-2 do begin
		temp1=strmid(temp, 0, strpos(temp, ","))
		mfheaderDescriptor.fabs[n].idxtype[l]=long(temp1)
		temp=strmid(temp, strpos(temp, ",")+1)
	endfor
	temp1=strmid(temp, 0, strpos(temp, ")"))
	mfheaderDescriptor.fabs[n].idxtype[l]=long(temp1)	
endfor

; read final parenthesis
readf, fp, temp

; read number of fabs stored in this set of multifab files
; this always appears to be the same as the number of fabs above,
; so I've coded this assuming it's redundant. Perhaps I will have to
; change this to accomodate more general file formats in the future.
readf, fp, temp

; now read the location of each fab
for n=0, nfab-1 do begin
	readf, fp, temp			; read line of text
	temp=strmid(temp, strpos(temp, ":")+2)		; skip junk
	mfheaderDescriptor.fabs[n].filename = $
		strmid(temp, 0, strpos(temp, " "))
						; portion before space
	temp=strmid(temp, strpos(temp, " ")+1)	; portion after space
	mfheaderDescriptor.fabs[n].offset = long(temp)
endfor

; read the min of each component in each fab

readf, fp, temp		; read a blank line
readf, fp, temp		; read the number of fabs and components again
for n=0, nfab-1 do begin
	readf, fp, temp		; read line to be parsed
	for m=0, ncomp-1 do begin
		temp1=strmid(temp, 0, strpos(temp, ","))
					; portion before comma
		mfheaderDescriptor.fabs[n].min[m]=double(temp1)
		temp=strmid(temp, strpos(temp, ",")+1)
					; skip to after next comma
	endfor
endfor

; repeat for max of each fab

readf, fp, temp		; read a blank line
readf, fp, temp		; read the number of fabs and components again
for n=0, nfab-1 do begin
	readf, fp, temp		; read line to be parsed
	for m=0, ncomp-1 do begin
		temp1=strmid(temp, 0, strpos(temp, ","))
					; portion before comma
		mfheaderDescriptor.fabs[n].max[m]=double(temp1)
		temp=strmid(temp, strpos(temp, ",")+1)
					; skip to after next comma
	endfor
endfor

; close file
free_lun, fp


return, mfheaderDescriptor

end