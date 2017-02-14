; uses a voxel, which is then rotated and projected onto a grid

function vox_to_image,vox,angles,cubic=cubic,imgsize=imgsize

if keyword_set(cubic) then begin
    cu = cubic
endif else begin
    cu = 0.0
endelse
sizevox = size(vox)
sizeangles = size(angles)


if sizevox[0] ne 3 then begin
    print,'Error in VOX_TO_IMAGE.  Voxel variable is unsuitable.'
endif

if (sizeangles[0] ne 1) or (sizeangles[1] ne 3) then begin
    print,'Error in VOX_TO_IMAGE.  Angles variable is unsuitable.'
endif
azimuth = angles[0]
altitude = angles[1]
pitch = angles[2]
stipvox = vox ; make local copy to manipulate

if azimuth ne 0.0 then begin
    for i=0,sizevox[3]-1 do stipvox[*,*,i] = $
      rot(reform(stipvox[*,*,i]),azimuth/!DTOR,/interp, cubic=cu)
endif
if altitude ne 0.0 then begin
    for i=0,sizevox[1]-1 do stipvox[i,*,*] = $
      rot(reform(stipvox[i,*,*]),altitude/!DTOR,/interp, cubic=cu)
endif
if pitch ne 0.0 then begin
    for i=0,sizevox[3]-1 do stipvox[*,*,i] = $
      rot(reform(stipvox[*,*,i]),pitch/!DTOR,/interp, cubic=cu)
endif

img = stipvox[*,*,0]
img[*] = 0.0

for i=0,sizevox[3]-1 do begin
    img = img + stipvox[*,*,i]
endfor

tvscl,img

return,img
end
