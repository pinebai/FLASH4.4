function slice_any, tree,angles,target=target, $
                     xrange=xrange, yrange=yrange, 	$
                     verbose=verbose, scalerange=scalerange,			$
                     maxlevel=maxlevel, log=log, $
                     subsample=subsample, $
                     imgdim=imgdim, $
                     pltsize=pltsize, rhocut=rhocut, smearsize=smearsize,	$
                     smeardist=smeardist, _extra = extra

; This program produces a slice of the data
; The slice plane is centered on target
; The plane is normal to Z, but rotated by 3 angles, azimuth, then
; altitude, then pitch
; xrange and yrange are the extent of the plot in the plane of the
; slice; transpose camera target
if not keyword_set(target) then target = [0.0,0.0,0.0]
azimuth = angles[0]
altitude = angles[1]
pitch = angles[2]

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


QSC = QR##(QP##(QA##QT)) ; quaternion for sim to camer
print,QSC
QCS = invert(QSC)

if not keyword_set(imgdim) then begin
    imgdim = [512,512]
endif
img = dblarr(imgdim[0],imgdim[1])


;for each fab begin
;    if M##fabcenter is close to plane then begin
;        if fab intersects plan then begin
;            ; optionally make pared-down list of cells close to plane
;            for each cell begin
;                if cell center close to plane then begin
;                    if cell intersects plane then begin
;                        s = set of indices of slice that might be in cell
;                        for each image[s] begin
;                            if image[s] is in cell, image[s] = cellvalue
;                        endfor
;                    endif
;                endif
;            endfor
;        endif
;    endif
;endfor
xlim = dblarr(2)
ylim = xlim
zlim = ylim
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
zlim = dblarr(2)
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

slice_node,tree.rootnode,QSC,img,[[xlim],[ylim],[zlim]],tree.gridspacing,amrlim
return,img
end

