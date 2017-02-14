function read_fabcomponent, fabnum, component, $
	plotDescriptor, mfheaderDescriptor

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


; this routine reads a single component of a fab. The component is
; specified as component, and the fab is specified as fabnum. The header
; for the set of multifabs in which this component is contained must
; already have been read inot mfheaderDescriptor.


; do error checking
if fabnum ge mfheaderDescriptor.nfab then begin
	print, 'Invalid fab number: ', fabnum
	print, 'Maximum fab number is ', mfheaderDescriptor.nfab-1
	return, -1
endif
if component ge mfheaderDescriptor.ncomp then begin
	print, 'Invalid component number: ', component
	print, 'Maximum component number is ', mfheaderDescriptor.ncomp-1
	return, -1
endif

; create an array to receive the output
if plotDescriptor.ndim eq 1 then $
	fabarr=dblarr(mfheaderDescriptor.fabs[fabnum].idxhi[0] - $
		mfheaderDescriptor.fabs[fabnum].idxlo[0]+1) $
else if plotDescriptor.ndim eq 2 then $
	fabarr=dblarr(mfheaderDescriptor.fabs[fabnum].idxhi[0] - $
		mfheaderDescriptor.fabs[fabnum].idxlo[0]+1, $
		mfheaderDescriptor.fabs[fabnum].idxhi[1] - $
		mfheaderDescriptor.fabs[fabnum].idxlo[1]+1) $
else if plotDescriptor.ndim eq 3 then $
	fabarr=dblarr(mfheaderDescriptor.fabs[fabnum].idxhi[0] - $
		mfheaderDescriptor.fabs[fabnum].idxlo[0]+1, $
		mfheaderDescriptor.fabs[fabnum].idxhi[1] - $
		mfheaderDescriptor.fabs[fabnum].idxlo[1]+1, $
		mfheaderDescriptor.fabs[fabnum].idxhi[2] - $
		mfheaderDescriptor.fabs[fabnum].idxlo[2]+1) $
else return, -1
fabsize=n_elements(fabarr)

; construct the name of the relevant multifab file
fname = plotDescriptor.name + "/" + mfheaderDescriptor.basename + $
	strmid(mfheaderDescriptor.fabs[fabnum].filename, $
	strpos(mfheaderDescriptor.fabs[fabnum].filename, "_"))

; open the file and point to the appropriate point
openr, fp, fname, /get_lun
point_lun, fp, mfheaderDescriptor.fabs[fabnum].offset

; read the header info 
temp=''
readf, fp, temp

; The first part of the header is information on how real numbers are
; stored in binary on the platform that generated this plot file.
; I'm not sure what all this stuff is, but the first number appears to
; be the number of bytes used to store each real number. Grab that.
temp=strmid(temp, strpos(temp, "((")+2)
temp1=strmid(temp, 0, strpos(temp, ","))
nbytes=long(temp1)

; The next thing we need is the byte ordering. From what I can tell,
; the object in the second set of parentheses is always
; (nbytes, (1 2 3 ... nbytes)) or
; (nbytes, (nbytes nbytes-1 nbytes-2 ... 1)).
; The first of these means little endian (same as the Sun and Intel),
; the second means big endian.

temp=strmid(temp, strpos(temp, "))")+2)
temp=strmid(temp, strpos(temp, "(")+1)
temp1=strmid(temp, 0, strpos(temp, ","))
nbytes2=long(temp1)
temp=strmid(temp, strpos(temp, "(")+1)
temp1=strmid(temp, 0, strpos(temp, " "))
firstnum=long(temp1)
if firstnum ne nbytes2 then bigendian=1 else bigendian=0

; This uses the non-standard routine is_ieee_big to determine the
; endian-ness of the local machine. A copy of the routine is included
; with amrlib, but I didn't write it.

if bigendian ne is_ieee_big() then do_swap_endian = 1 else do_swap_endian = 0

; Next are the index limits on this fab. Since these should already be
; recorded in the multifab header descriptor, it should be safe to skip
; them too.

; At this point, the file pointer should be at the start of the binary
; data for the fab. The file is arranged with all of the first
; component, then all of the second component, etc. Now, skip to the
; appropriate component.
point_lun, -fp, posn		; get current position
point_lun, fp, posn + component * nbytes * fabsize
				; jump to start of component we want

; Now read the component
readu, fp, fabarr

; Swap endianness if needed
if do_swap_endian then fabarr = swap_endian(fabarr)

; Close file and return
free_lun, fp
return, fabarr

end