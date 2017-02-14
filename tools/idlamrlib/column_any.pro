function column_any, tree,angles,target=target,anisotropic=anositropic, $
                     xrange=xrange, yrange=yrange, zrange=zrange, nocolorbar=nocolorbar,	$
                     verbose=verbose, scalerange=scalerange,			$
                     maxlevel=maxlevel, log=log, snaptogrid=snaptogrid,	$
                     subsample=subsample, $
                     imgdim=imgdim, $
                     voxdim=voxdim, $
                     pltsize=pltsize, rhocut=rhocut, smearsize=smearsize,	$
                     smeardist=smeardist, blur=blur, _extra = extra

; This program produces a column density plot of an AMR-oct-tree object as
; viewed from an arbitrary angle.  Otherwise, it behaves like
; column_amr
; Make the tree object by using amr_tree
; The argument angles is a fltarr(3) expressing the orientation of the
; AMR object in the 3D window.
; angles = [altitude,azimuth,pitch]
; Example: [0,0,0] is like looking at the North Pole
;          [Pi/2,0,0] is like looking at West Africa
;          [Pi/2,0,Pi/2] is like looking at West Africa, rotated so
;          Europe it to the right.
;          [2*Pi/3,2*Pi/3,0] is like looking at Australia
;          [0,0,Pi/2] is like rotating clockwise about the north pole
; target is the location in the AMR object that the camera is centered
; on.  It coordinates are the normal XYZ.  If omitted it defaults to
; the center of the domain.  target is given in AMR coords
; xrange, yrange,zrange determine what parts of the amr object get
; rendered
; the x,y,z values are in camera space, not amr space
; (*((*(tree.rootnode.children[i])).children[*])).value
print,systime()
azimuth = angles[0]
altitude = angles[1]
pitch = angles[2]
LocValue = {TypeLocValue, x:fix(0), y:fix(0), value:double(0.0)}

ViewParam = {$
              amrfilename:'', $
              target:fltarr(3),$
              altitude:0.0,azimuth:0.0,pitch:0.0,$
              xrange:fltarr(2),yrange:fltarr(2),zrange:fltarr(2),$
              scalerange:fltarr(2),$
              maxlevel:-1}

if not keyword_set(target) then target = (tree.boxmin + tree.boxmax) / 2.0


; the identity quaternion
QI = [[1,0,0,0],$
      [0,1,0,0],$
      [0,0,1,0],$
      [0,0,0,1]]

; transpose camera target
QT = [[1,0,0,-target[0]],$
      [0,1,0,-target[1]],$
      [0,0,1,-target[2]],$
      [0,0,0,1]]

; change camera azimuth angle
QA = [[ cos(azimuth),sin(azimuth),0,0],$
      [-sin(azimuth),cos(azimuth),0,0],$
      [            0,           0,1,0],$
      [            0,           0,0,1]]

; change camera polar angle
QP = [[1,            0,             0,0],$
      [0, cos(altitude),sin(altitude),0],$
      [0,-sin(altitude),cos(altitude),0],$
      [0,             0,            0,1]]

; roll the camera angle
QR = [[ cos(pitch),sin(pitch),0,0],$
      [-sin(pitch),cos(pitch),0,0],$
      [          0,         0,1,0],$
      [          0,         0,0,1]]


QSC = QR##(QP##(QA##QT)) ; quaternion for sim to camera

QCS = invert(QSC)        ; quaternion for camera to sim

Q = QSC

; cameraspace = Q ## simulation space

if not keyword_set(outputmethod) then outputmethod = 'x'


outD = '' ; flag for set_plot,outD
outS = '' ; suffix for filename, if used
scalepix = 0 ; does plotting device have scalable pixels
if outputmethod eq 'x' then begin
    outM = 'x'
    outD = 'X'
    altout = 'X'
    outS = ''
    scalepix = 0
endif else if outputmethod eq 'mac' then begin
    outM = 'mac'
    outD = 'MAC'
    altout = 'MAC'
    outS = ''
    scalepix = 0
endif else if outputmethod eq 'win' then begin
    outM = 'win'
    outD = 'WIN'
    altout = 'WIN'
    outS = ''
    scalepix = 0
endif else if outputmethod eq 'z' then begin
    outM = 'z'
    outD = 'Z'
    altout = 'X'
    outS = ''
    scalepix = 0
endif else if outputmethod eq 'ps' then begin
    outM = 'ps'
    outD = 'PS'
    altout = 'X'
    outS = '.ps'
    scalepix = 1
endif else if outputmethod eq 'cgm' then begin
    outM = 'cgm'
    outD = 'CGM'
    altout = 'X'
    outS = '.cgm'
    scalepix = 1
endif else if outputmethod eq 'ppm' then begin
    outM = 'ppm'
    outD = 'Z'
    altout = 'X'
    outS = '.ppm'
    scalepix = 0
endif else if outputmethod eq 'jpeg' then begin
    outM = 'jpeg'
    outD = 'Z'
    altout = 'X'
    outS = '.jpeg'
    scalepix = 0
endif


xlim = dblarr(2)
ylim = dblarr(2)
zlim = dblarr(2)
if keyword_set(xrange) then begin
    xlim=xrange
endif else begin
    xlim[0] = min( [ $
                     (QSC##[tree.boxmin[0],0,0,1])[0], $
                     (QSC##[0,tree.boxmin[1],0,1])[0], $
                     (QSC##[0,0,tree.boxmin[2],1])[0], $
                     (QSC##[tree.boxmax[0],0,0,1])[0], $
                     (QSC##[0,tree.boxmax[1],0,1])[0], $
                     (QSC##[0,0,tree.boxmax[2],1])[0]  $
                   ] )
    xlim[1] = max( [ $
                     (QSC##[tree.boxmin[0],0,0,1])[0], $
                     (QSC##[0,tree.boxmin[1],0,1])[0], $
                     (QSC##[0,0,tree.boxmin[2],1])[0], $
                     (QSC##[tree.boxmax[0],0,0,1])[0], $
                     (QSC##[0,tree.boxmax[1],0,1])[0], $
                     (QSC##[0,0,tree.boxmax[2],1])[0]  $
                   ] )
endelse

if keyword_set(yrange) then begin
    ylim=yrange
endif else begin
    ylim[0] = min( [ $
                     (QSC##[tree.boxmin[0],0,0,1])[1], $
                     (QSC##[0,tree.boxmin[1],0,1])[1], $
                     (QSC##[0,0,tree.boxmin[2],1])[1], $
                     (QSC##[tree.boxmax[0],0,0,1])[1], $
                     (QSC##[0,tree.boxmax[1],0,1])[1], $
                     (QSC##[0,0,tree.boxmax[2],1])[1]  $
                   ] )
    ylim[1] = max( [ $
                     (QSC##[tree.boxmin[0],0,0,1])[1], $
                     (QSC##[0,tree.boxmin[1],0,1])[1], $
                     (QSC##[0,0,tree.boxmin[2],1])[1], $
                     (QSC##[tree.boxmax[0],0,0,1])[1], $
                     (QSC##[0,tree.boxmax[1],0,1])[1], $
                     (QSC##[0,0,tree.boxmax[2],1])[1]  $
                   ] )
endelse

if keyword_set(zrange) then begin
    zlim=zrange
endif else begin
    zlim[0] = min( [ $
                     (QSC##[tree.boxmin[0],0,0,1])[2], $
                     (QSC##[0,tree.boxmin[1],0,1])[2], $
                     (QSC##[0,0,tree.boxmin[2],1])[2], $
                     (QSC##[tree.boxmax[0],0,0,1])[2], $
                     (QSC##[0,tree.boxmax[1],0,1])[2], $
                     (QSC##[0,0,tree.boxmax[2],1])[2]  $
                   ] )
    zlim[1] = max( [ $
                     (QSC##[tree.boxmin[0],0,0,1])[2], $
                     (QSC##[0,tree.boxmin[1],0,1])[2], $
                     (QSC##[0,0,tree.boxmin[2],1])[2], $
                     (QSC##[tree.boxmax[0],0,0,1])[2], $
                     (QSC##[0,tree.boxmax[1],0,1])[2], $
                     (QSC##[0,0,tree.boxmax[2],1])[2]  $
                   ] )
endelse

amrlimlo = dblarr(3)
amrlimhi = dblarr(3)

for i=0,2 do begin
    amrlimlo[i] = min( [ $
                         (QCS##[xlim[0],0,0,1])[i], $
                         (QCS##[0,ylim[0],0,1])[i], $
                         (QCS##[0,0,zlim[0],1])[i], $
                         (QCS##[xlim[1],0,0,1])[i], $
                         (QCS##[0,ylim[1],0,1])[i], $
                         (QCS##[0,0,zlim[1],1])[i]  $
                       ] )
    amrlimhi[i] = max( [ $
                         (QCS##[xlim[0],0,0,1])[i], $
                         (QCS##[0,ylim[0],0,1])[i], $
                         (QCS##[0,0,zlim[0],1])[i], $
                         (QCS##[xlim[1],0,0,1])[i], $
                         (QCS##[0,ylim[1],0,1])[i], $
                         (QCS##[0,0,zlim[1],1])[i]  $
                       ] )
endfor
amrlim = [[amrlimlo],[amrlimhi]]
print,amrlim


;if not keyword_set(xrange) then begin
;    ;xlim=[amr.boxmin[dim1], amr.boxmax[dim1]] 
;    xlim = [ (Q##[tree.boxmin,1.0])[0], (Q##[tree.boxmax,1.0])[0] ]
;    if xlim[0] gt xlim[1] then xlim = reverse(xlim)
;endif else xlim=xrange
;if not keyword_set(yrange) then begin
;    ylim= [ (Q##[tree.boxmin,1.0])[1], (Q##[tree.boxmax,1.0])[1] ]
;    if ylim[0] gt ylim[1] then ylim = reverse(ylim)
;endif else ylim=yrange
;if not keyword_set(zrange) then begin
;    zlim= [ (Q##[tree.boxmin,1.0])[2], (Q##[tree.boxmax,1.0])[2] ]
;    if zlim[0] gt zlim[1] then zlim = reverse(zlim)
;endif else zlim = zrange










; decide on number of pixels
; based on # pixels,zoom,ref ratios, pick finest level to render
if not scalepix then begin
    ; we do not have scalable pixels
    ; get max and min for image box in pixels
    mindx = tree.gridspacing[0,tree.maxlevel-1]
    imgdisplo = convert_coord(xlim[0]/mindx, ylim[0]/mindx, /data, /to_device)
    imgdisphi = convert_coord(xlim[1]/mindx, ylim[1]/mindx, /data, /to_device)
    imgdisplo=round(imgdisplo[0:1])
    imgdisphi=round(imgdisphi[0:1])
    if (imgdisphi[0] - imgdisplo[0]) lt 2 then begin
        print, 'Error: requested x axis range is too small.'
        print,qsc
        ;blargh
    endif
    if (imgdisphi[1] - imgdisplo[1]) lt 2 then begin
        print, 'Error: requested y axis range is too small.'
        print,qsc
        ;blargh
    endif    
endif
if not keyword_set(imgdim) then begin
    imgdim = [512,512,1]
endif
if keyword_set(subsample) then begin
    ; make an image subsample times as
    ; large and then average back to smaller image
    imgdim = imgdim * subsample
endif else begin
    subsample = 1.0
endelse
imgdimx = imgdim[0]
imgdimy = imgdim[1]
imgdimz = imgdim[2] ; not applicable for 2-d images, but should still be set
img = dblarr(imgdimx,imgdimy,imgdimz)
img[*] = 0.0

if not keyword_set(maxlevel) then begin
    if scalepix then begin
        ml = tree.maxlevel
    endif else begin
        ml = min([tree.maxlevel, fix(alog10(imgdimx/32.0)/alog10(4.0) + 0.5)])
    endelse
    maxlevel=ml
endif else ml=maxlevel

if scalepix and not keyword_set(maxlevel) then maxlevel=finestlevel


if not keyword_set(anisotropic) then isotropic=1 else isotropic=0
if keyword_set(maxlevel) then begin
    ml=maxlevel
endif else begin
    ; decide the max level that can be resolved
    mdx = (xlim[1] - xlim[0]) / imgdimx
    mdy = (ylim[1] - ylim[0]) / imgdimy
    mlx = max(where(tree.gridspacing[0,*] gt mdx)) 
    mly = max(where(tree.gridspacing[0,*] gt mdy)) 
    if mlx gt mly then begin
        ml = mlx 
    endif else ml = mly
    if ml gt tree.maxlevel then begin
        ml=tree.maxlevel
    endif
endelse

; go through the tree 
;  if a whole fab is smaller then mdx, subsample the fab (subsample)^3
;  times and add to the correct array point(s)
;  if a fab is bigger than mdx, for each cell in it, subsample is and
;  add the cell parts to array points
;  if a node is smaller than mdx, do subsample and add to array, and
;  do not recurse
;img = img + project_node($
;                              tree.rootnode,$
;                              QSC,$
;                              [xlim,ylim,zlim],$
;                              img, $
;                              maxlevel, $
;                              subsample $
;                            )

stipsize = intarr(ml+1, 2)
stipptrs = ptrarr(ml+1)
idealstipsize = dblarr(ml+1)
;for l=0,ml do begin
;endfor
print,maxlevel,ml
for l=0,ml do begin
    ; tempsize = maximum extent of a cell in any one dimension
    ; for cubes, it is sqrt(3)*cellsize
  tempsize = 2.0* sqrt( $
                        (tree.gridspacing[0,l])^2 + $
                        (tree.gridspacing[1,l])^2 + $
                        (tree.gridspacing[2,l])^2 )
  stipsize[l,*] = [ $
                    tempsize / double((xlim[1]-xlim[0]) / double(imgdimx)), $
                    tempsize / double((ylim[1]-ylim[0]) / double(imgdimy)) ]
  idealstipsize[l] = max([tempsize / double((xlim[1]-xlim[0]) / double(imgdimx)),$
                       tempsize / double((ylim[1]-ylim[0]) / double(imgdimy))])
  ;stipsize = stipsize * 2.0
  if stipsize[l,0] lt 1 then stipsize[l,0]=1
  if stipsize[l,1] lt 1 then stipsize[l,1]=1

  cu = 0.0 ; cubic-interpolation parameter 
  if (stipsize[l,0] gt 1) or (stipsize[l,1] gt 1) then begin
      if (double(stipsize[l,0])*double(stipsize[l,1])) ge double(2.0)^24 then begin
          ; voxel so huge just the mask might blow memory
          ; just make something large and blury
          mss = 4096
          stipptrs[l] = ptr_new(dblarr(mss,mss))
          (*(stipptrs[l]))[*] = double(1.0)/(double(2.0)^24)
          print,'Warning: Could not make '+string(stipsize[l,0])+' by '+string(stipsize[l,1])+' mask.  Using bogus 4096^2 mask.'
      endif else if (stipsize[l,0] gt 100) or (stipsize[l,1] gt 100) then begin
          ; don't attempt huge voxels
          ; scale down the problem, then scale it back up
          mss = 100 
          stipvox = dblarr(100,100,100) 
          cellen = tree.gridspacing[*,l]
          cellen = cellen / sqrt((cellen[0])^2+ (cellen[1])^2 + (cellen[2])^2  )
          stipdx = 100*cellen[0]
          stipdy = 100*cellen[1] 
          stipdz = 100*cellen[2]
          msss = fix( [(mss-stipdx)/2.0,(mss-stipdy)/2.0, (mss-stipdz)/2.0 ])
                                ; stipdx = xcellsize / xpixsize
          stipvox[ msss[0] : msss[0]+stipdx-1, $
                   msss[0] : msss[1]+stipdy-1, $
                   msss[0] : msss[2]+stipdz-1 ] = 1.0
          if azimuth ne 0.0 then begin
              for i=0,mss-1 do stipvox[*,*,i] = $
                rot(reform(stipvox[*,*,i]),azimuth/!DTOR,missing=0,/interp, cubic=cu)
          endif
          if altitude ne 0.0 then begin
              for i=0,mss-1 do stipvox[i,*,*] = $
                rot(reform(stipvox[i,*,*]),altitude/!DTOR,missing=0,/interp, cubic=cu)
          endif
          if pitch ne 0.0 then begin
              for i=0,mss-1 do stipvox[*,*,i] = $
                rot(reform(stipvox[*,*,i]),pitch/!DTOR,missing=0,/interp, cubic=cu)
          endif
          ; correct from any neg values introduced from interpolation
          wh = where(stipvox lt 0.0)
          if wh[0] ne -1 then stipvox[wh] = 0.0

          stipptrs[l] = ptr_new(dblarr(stipsize[l,0],stipsize[l,1]))
          *(stipptrs[l]) = 0.0
          ;for i=0,mss-1 do *(stipptrs[l]) = *(stipptrs[l]) + congrid(reform(stipvox[*,*,i]),stipsize[l,0],stipsize[l,1], 1,/interp)
          tempstip = dblarr(100,100)
          for i=0,99 do tempstip = tempstip + reform(stipvox[*,*,i])
          *(stipptrs[l]) = congrid(tempstip,stipsize[l,0],stipsize[l,1],/interp, cubic=cu)
          *(stipptrs[l]) = *(stipptrs[l]) / total(*(stipptrs[l]) ) 
          ;bork
      endif else begin
          mss = max([stipsize[l,0], stipsize[l,1]])
          ;mss = sqrt(stipsize[l,0]^2 + stipsize[l,1]^2 + (stipsize[l,0]+stipsize[l,1])^2/4.0)
          ;mss = round(mss)
          if (mss mod 2 ) eq 1 then mss = mss+1 ; for better results with ROT(...)
          stipvox = dblarr(mss,mss,mss) 
          stipdx = tree.gridspacing[2,l] / double(((xlim[1]-xlim[0])) / double(imgdimx))
          stipdy = tree.gridspacing[2,l] / double(((ylim[1]-ylim[0])) / double(imgdimy))
          print,'StipDXDY',stipdx,stipdy
          if stipdx ne stipdy then print,bork
          stipdz = max([stipdx,stipdy])

          if stipdx gt 1.0 then begin
              stedx = stipdx - fix(stipdx)
              stipdx = fix(stipdx)
          endif else begin
              stedx = 0.0
              stipdx = 1.0
          endelse

          if stipdy gt 1.0 then begin
              stedy = stipdy - fix(stipdy)
              stipdy = fix(stipdy)
          endif else begin
              stedy = 0.0
              stipdy = 1.0
          endelse

          if stipdz gt 1.0 then begin
              stedz = stipdz - fix(stipdz)
              stipdz = fix(stipdz)
          endif else begin
              stedz = 0.0
              stipdz = 1.0
          endelse
          ;stedx = 0 & stedy = 0 & stedz = 0
          
          msss = fix( [(mss-stipdx)/2.0,(mss-stipdy)/2.0, (mss-stipdz)/2.0 ])
                                ; stipdx = xcellsize / xpixsize
          print,stedx,stedy,stedz
          ; fill center
          stipvox[ msss[0] : msss[0]+stipdx-1, $
                   msss[1] : msss[1]+stipdy-1, $
                   msss[2] : msss[2]+stipdz-1 ] = 1.0

          ; fill edge planes
          if stedx gt 0.0 then begin
              stipvox[ msss[0]-1, $ ; fill low-x plane
                       msss[1] : msss[1]+stipdy-1, $
                       msss[2] : msss[2]+stipdz-1 ] = stedx/2.0
              stipvox[ msss[0]+stipdx, $ ; fill high-x plane
                       msss[1] : msss[1]+stipdy-1, $
                       msss[2] : msss[2]+stipdz-1 ] = stedx/2.0
          endif
          if stedy gt 0.0 then begin
              stipvox[ msss[0] : msss[0]+stipdx-1, $
                       msss[1]-1, $ ; fill low-y plane
                       msss[2] : msss[2]+stipdz-1 ] = stedy/2.0
              stipvox[ msss[0] : msss[0]+stipdx-1, $
                       msss[1]+stipdy, $ ; fill high-y plane
                       msss[0] : msss[2]+stipdz-1 ] = stedy/2.0
          endif
          if stedz gt 0.0 then begin
              stipvox[ msss[0] : msss[0]+stipdx-1, $
                       msss[1] : msss[1]+stipdy-1, $ ; fill low-z plane
                       msss[2]-1]                   = stedz/2.0
              stipvox[ msss[0] : msss[0]+stipdx-1, $
                       msss[1] : msss[1]+stipdy-1, $ ; fill high-z plane
                       msss[2]+stipdx]              = stedz/2.0
          endif

          ; fill edge lines
          if stedx+stedy gt 0.0 then begin
              stipvox[msss[0]-1, msss[1]-1, msss[0] : msss[2]+stipdz-1] = stedx*stedy*0.25 ; lo-x lo-y
              stipvox[msss[0]-1, msss[1]+stipdy, msss[0] : msss[2]+stipdz-1] = stedx*stedy*0.25 ; lo-x hi-y
              stipvox[ msss[0]+stipdx, msss[1]-1, msss[0] : msss[2]+stipdz-1] = stedx*stedy*0.25 ; hi-x lo-y
              stipvox[ msss[0]+stipdx, msss[1]+stipdy, msss[0] : msss[2]+stipdz-1] = stedx*stedy*0.25 ; hi-x hi-y
          endif
          if stedy+stedz gt 0.0 then begin
              stipvox[msss[0] : msss[0]+stipdx-1, msss[1]-1, msss[2]-1] = stedy*stedz*0.25 ; lo-y lo-z
              stipvox[msss[0] : msss[0]+stipdx-1, msss[1]-1, msss[2]+stipdz] = stedy*stedz*0.25 ;lo-y hi-z              
              stipvox[msss[0] : msss[0]+stipdx-1, msss[1]+stipdy, msss[2]-1] = stedy*stedz*0.25 ; hi-y lo-z
              stipvox[msss[0] : msss[0]+stipdx-1, msss[1]+stipdy, msss[2]+stipdz] = stedy*stedz*0.25 ;hi-y hi-z
          endif
          if stedz+stedx gt 0.0 then begin          
              stipvox[msss[0]-1, msss[1] : msss[1]+stipdy-1,msss[2]-1] = stedx*stedz*0.25 ; lo-x lo-z
              stipvox[msss[0]-1, msss[1] : msss[1]+stipdy-1,msss[2]+stipdz] = stedx*stedz*0.25 ; lo-x hi-z          
              stipvox[ msss[0]+stipdx, msss[1] : msss[1]+stipdy-1,msss[2]-1] = stedx*stedz*0.25 ; hi-x lo-z
              stipvox[ msss[0]+stipdx, msss[1] : msss[1]+stipdy-1,msss[2]+stipdz] = stedx*stedz*0.25 ; hi-x hi-z
          endif
          if stedx+stedy+stedz gt 0.0 then begin                            
                                ; fill corners
              stipvox[msss[0]-1, msss[1]-1, msss[2]-1] = stedx*stedy*stedz*0.125
              stipvox[msss[0]-1, msss[1]-1, msss[2]+stipdz] = stedx*stedy*stedz*0.125
              stipvox[msss[0]-1, msss[1]+stipdy, msss[2]-1] = stedx*stedy*stedz*0.125
              stipvox[msss[0]-1, msss[1]+stipdy, msss[2]+stipdz] = stedx*stedy*stedz*0.125
              stipvox[msss[0]+stipdx, msss[1]-1, msss[2]-1] = stedx*stedy*stedz*0.125
              stipvox[msss[0]+stipdx, msss[1]-1, msss[2]+stipdz] = stedx*stedy*stedz*0.125
              stipvox[msss[0]+stipdx, msss[1]+stipdy, msss[2]-1] = stedx*stedy*stedz*0.125
              stipvox[msss[0]+stipdx, msss[1]+stipdy, msss[2]+stipdz] = stedx*stedy*stedz*0.125
          endif

          ;print,total(stipvox)
          stipptrs[l] = ptr_new(dblarr(mss,mss))
          if azimuth ne 0.0 then begin
              for i=0,mss-1 do stipvox[*,*,i] = $
                rot(reform(stipvox[*,*,i]),azimuth/!DTOR,/interp, cubic=cu)
          endif
          if altitude ne 0.0 then begin
              for i=0,mss-1 do stipvox[i,*,*] = $
                rot(reform(stipvox[i,*,*]),altitude/!DTOR,/interp, cubic=cu)
          endif
          if pitch ne 0.0 then begin
              for i=0,mss-1 do stipvox[*,*,i] = $
                rot(reform(stipvox[*,*,i]),pitch/!DTOR,/interp, cubic=cu)
          endif
          ; correct for any neg values introduced by interpolation
          wh = where(stipvox lt 0.0)
          if (wh[0] ne -1) then stipvox[wh] = 0.0

          for i=0,mss-1 do *(stipptrs[l]) = *(stipptrs[l]) + reform(stipvox[*,*,i])
          if ( total(*(stipptrs[l])) gt 0.0) then begin
              *(stipptrs[l]) = *(stipptrs[l]) / total(*(stipptrs[l]) ) 
          endif else begin
              print,'Error, stipptrs[l] has zero values!'
              bork
          endelse
      endelse
  endif else begin
      stipptrs[l] = ptr_new(dblarr(1,1)+1.0)
  endelse
  sx = (size(*(stipptrs[l])))[1]
  sy = (size(*(stipptrs[l])))[2]
  ;if sx gt 3 then begin
  ;    kern = sx/4 - dist(sx/2,1)
  ;    ;exp( -1*(sx/8/sigma)^2) = 0.5
  ;    ;(sx/8/sigma)^2 = 0.7
  ;    ;sx/8/sigma = 0.832555
  ;    ;sigma = sx/8/0.83255
  ;    kern = exp(-1*(kern/2.0/0.83255)^2)
  ;    kern = kern#kern
  ;    *(stipptrs[l]) = convol(*(stipptrs[l]),kern)
  ;endif

  if keyword_set(blur) then begin
      ;print,'BLUR!!!!!!!!!!!!!!'
      if ((sx ge 12) and (sy ge 12)) then begin
          *(stipptrs[l]) = smooth(*(stipptrs[l]),sx/4)
      endif else if ((sx ge 3) and (sy ge 3)) then *(stipptrs[l]) = smooth(*(stipptrs[l]),3)
  endif
  a = min( where((*(stipptrs[l]))) / sy )-1
  b = max( where((*(stipptrs[l]))) / sy )+1
  c = min( where((*(stipptrs[l]))) mod sy)-1
  d = max( where((*(stipptrs[l]))) mod sy)+1
  if a lt 0 then a = 0
  if b ge stipsize[l,1] then b = stipsize[l,1]-1
  if c lt 0 then c = 0
  if d ge stipsize[l,0] then d = stipsize[l,0]-1
  ;stipsize[l,*] = [ d-c+1, b-a+1  ] 
  if ((b-a) lt 4096 and (d-c lt 4096) ) and (not ((b-a gt 0) xor (d-c gt 0))) then begin
      ;trim the stipple mask down if it is not 1-by-X or X-by-1 in size
      newstip = (*(stipptrs[l]))[c:d,a:b]
      ptr_free,stipptrs[l]
      stipptrs[l] = ptr_new(newstip)      
      stipsize[l,*] = [ d-c+1, b-a+1  ] 
  endif else begin
      stipsize[l,*] = [mss,mss]
  endelse
  if total(*(stipptrs[l])) ne double(1.0) then begin
  ;    print,'Error, total is now '+string(total(*(stipptrs[l])))
  ;    print,bork
      *(stipptrs[l]) = *(stipptrs[l]) / total(*(stipptrs[l]) )
  endif
  tv_colorized,bytscl(congrid(*(stipptrs[l]),imgdimx,imgdimy,1)) & print,systime()
  ; trim stipptrs[l] to be as small as possible
  ;stipsize[l,0],stipsize[l,1]
  ;(*(stipptrs[0]))[c:d,a:b]
  
  tv_colorized,bytscl(*(stipptrs[l]) )
  ;bork
  ;wait,1
endfor
print,'Stipsize = '
print,stipsize

print,'Calling PROJECT_NODE '
print,'Node = ',tree.rootnode
print,'Q = ',QSC
print,'Size(img) = ',size(img)
print,'Imglim = ',[[xlim],[ylim],[zlim]]
print,'Maxlevel = ',ml
print,'Subsample = ',subsample
print,'Gridspacing = ',tree.gridspacing
project_node,tree.rootnode,QSC,img,[[xlim],[ylim],[zlim]],ml,subsample,tree.gridspacing,stipptrs,amrlim
if keyword_set(subsample) then begin
    imgdim = imgdim / subsample
    img = congrid(img,imgdim[0],imgdim[1],imgdim[2],/interp)
endif
;tvscl,alog10(img)
imt = img
whi = where(imt gt 0.0)
if total(whi) ne -1 then begin
    imm = min(imt[whi])
endif else imm = 0.0
imw = where(imt lt imm)
if total(imw) ne -1 then begin
    imt[imw] = imm
endif
imt = alog10(imt)
imt = bytscl(imt)
tv_colorized,imt
;loadct,13
;c = intarr(256,3)
;tvlct,c,/get
;imc = intarr(imgdimx,imgdimy,3)
;for i=0,2 do imc[*,*,i] = c[imt[*,*],i]
;tv,imc,true=3

print,'Image complete!'
help,img
print,'Stipsize = '+string(stipsize)
print,systime()
;blargh
if not keyword_set(voxdim) then begin
return,img
endif else begin
    data = {img:img,vox:voxel}
    return,data
endelse

end ; COLUMN_ANY.PRO

