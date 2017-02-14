pro project_node,node,QSC,img,imglim,maxlevel,subsample,gridspacing,stipptrs,amrlim
error=0

;if node.level gt maxlevel then begin
;    nodedim = (node.xhi-node.xlo)
;    nodeden = node.value / (double(nodedim[0])*double(nodedim[1])*double(nodedim[2]))
;    print,'Projecting node as cube with density = ',nodeden, nodedim
;    project_cube,img, imglim, QSC, nodeden, (node.xlo+node.xhi)/2.0,
;    (node.xhi-node.xlo), subsample

ximgsize = (size(img))[1]
yimgsize = (size(img))[2]
if (size(img))[0] eq 3 then begin
    zimgsize = (size(img))[3]
    usevox = 1
endif else begin
    zimgsize = 1
    usevox = 0
endelse
xcellsize = double((imglim[1,0]-imglim[0,0]))/ double(ximgsize)
ycellsize = double((imglim[1,1]-imglim[0,1]))/ double(yimgsize)
zcellsize = double((imglim[1,2]-imglim[0,2]))/ double(zimgsize)
QSCNR = QSC
QSCNR[3,0:2] = 0.0 ; set the transpose elements to zero

if 0 then begin
endif else begin

    for chin=0,7 do begin
        if ptr_valid(node.children[chin]) then begin
            ncc = *(node.children[chin])
            if n_tags(ncc,/length) eq n_tags({TypeNode},/length) then begin
                if 1 then begin ; (total( node.xhi lt amrlim[*,0]) + total(node.xlo gt amrlim[*,1])) eq 0 then begin
                    project_node,ncc,QSC,img,imglim,maxlevel,subsample,gridspacing,stipptrs,amrlim
                endif else begin
                                ; node falls outside camera space
                                ; projected onto amr domain
                endelse
            endif
        endif
    endfor

    for chi=0,7 do begin
        if ptr_valid(node.children[chi]) then begin
            ncc = *(node.children[chi])
            ;if n_tags(ncc) eq n_tags({TypeNode}) then begin
            ;    ; ncc is a node
            ;    ; PROJECT_NODE,ncc,QSC,img,imglim,maxlevel,subsample
            ;    ;ptr_free,borkborkbork
            ;    
            ;    ; check that child intersects imglim
            ;    if (total( node.xhi lt amrlim[*,0]) + total(node.xlo gt amrlim[*,1])) eq 0 then begin
            ;        project_node,ncc,QSC,img,imglim,maxlevel,subsample,gridspacing,stipptrs,amrlim
            ;    endif else begin
            ;      ; node falls outside camera space
            ;      ; projected onto amr domain
            ;    endelse
            if n_tags(ncc,/length) eq n_tags({TypeFab},/length) then begin
                ; ncc is a Fab
                ; if entire fab is small, do PROJECT_CUBE on it
                ; else PROJECT_CUBE on each cell
                tempsize = double(gridspacing[*,ncc.level])
                ir = (ncc.indhi - ncc.indlo) + [1,1,1]
                ir20 = ir[0]/2
                ir21 = ir[1]/2
                ir22 = ir[2]/2
                if (ncc.level eq -1) or (ncc.level gt maxlevel)  then begin
                    ; dummy statement to optionally turn off some levels
                endif else begin
                    ;if 0 gt 1 then begin
                    data = *(ncc.data)
                    VolOverArea = tempsize[0]*tempsize[1]*tempsize[2] / xcellsize / ycellsize
                    Vol1OverVol2 = tempsize[0]*tempsize[1]*tempsize[2] / xcellsize/ycellsize/zcellsize
                    if Vol1OverVol2 lt ((double(64.0)^(double(-4.0)))) then begin
                        print,'Too small!'
                    endif
                    
                    ; center of the fab, in camera space
                    fabcenterc = (QSC##[ncc.center,1.0])[0:2]
                    fabdxc = (QSCNR##[tempsize[0],0,0,1.0])[0:2]
                    fabdyc = (QSCNR##[0,tempsize[1],0,1.0])[0:2]
                    fabdzc = (QSCNR##[0,0,tempsize[2],1.0])[0:2]

                    ; code to check if fab intersects window
                    ; if not, rig outer for loop index so it aborts
                    fabmin = (fabcenterc - ir20*fabdxc - ir21*fabdyc - ir22*fabdzc)
                    fabmax = (fabcenterc + ir20*fabdxc + ir21*fabdyc + ir22*fabdzc)
                    simflo = ncc.center - [ir20*tempsize[0] , ir21*tempsize[1] , ir22*tempsize[2]]
                    simfhi = ncc.center + [ir20*tempsize[0] , ir21*tempsize[1] , ir22*tempsize[2]]
                    
                    intersects = 0
                    ;for inti = 0,2 do begin
                    ;    intersects = intersects + ((fabmin[inti] lt imglim[1,inti]) and (fabmax[inti] gt imglim[0,inti]))
                    ;endfor
                    for inti=0,2 do begin
                        intersects = intersects + ((simflo[inti] lt amrlim[inti,1]) and (simfhi[inti] gt amrlim[inti,0]))
                    endfor
                    if  intersects ne 3 then begin
                        fistart = ir[0]
                        ;print,'Not rendering fab centered at ',string(fabcenterc),intersects
                    endif else begin
                        ;print,'Rendering FAB ',intersects,ncc.indlo,ncc.indhi,ncc.value
                        fistart = 0
                    endelse
                    ;fistart=0

                    l = ncc.level
                    if l ge (size(stipptrs))[1] then print,'Stipsize too small uh oh',bork
                    stipsize = (size(*(stipptrs[l])))
                    if stipsize[0] eq 1 then begin
                        stipsize = [1,stipsize[1]]
                    endif else begin
                        stipsize = stipsize[1:2]
                    endelse

                    ; for randomized jitter to aleviate aliasing
                    sj = (tempsize[0]/double(xcellsize))
                    ; picking the amount by which to
                    ; jitter a cell is an art, not a science
                    if sj gt 0.0 then begin
                        ;sj = sj*sj / (sj + 0.0) / 8.0
                        ;sj = sj * (0.50 - 0.50/(1+(sj+0.0)^(0.50)))*0.5
                        sj = 0.0
                        ;sj = sj * 0.25
                    endif else begin
                        sj = 0.0
                    endelse
                    common seed,s 

                    ;fcpixx = ((QSC##[ncc.center,1.0])[0] -imglim[0,0])/xcellsize
                    ;fcpixy = ((QSC##[ncc.center,1.0])[1] -imglim[0,1])/ycellsize
                    ;cpixdx = (QSC##[0,tempsize[1],0,1.0])[0:1] / [xcellsize,ycellsize]
                    ;cpixdy = (QSC##[0,tempsize[1],0,1.0])[0:1] / [xcellsize,ycellsize]
                    ;cpixdz = (QSC##[0,0,tempsize[2],1.0])[0:2] / [xcellsize,ycellsize]
                    ;indbase = [fcpixx,fcpixy]
                    ;indplusx = (0.0-ir21 + 0.5 -1.0)*cpixdx
                    
                    cellcentcam = fabcenterc

                    ; init x-loop
                    relxcam = (0.0 - ir20 + 0.5 -1.0)*fabdxc
                    inci=1
                    incj = inci
                    inck = incj
                    for fi=fistart,ir[0]-1,inci do begin
                        ;relxind = fi - ir20 + 0.5
                        ;indplusx = indplusx + cpixdx
                        ;indplusxy 

                        

                        ; increment relxcam
                        relxcam = relxcam + fabdxc
                        ; update centplusx to be used in y-loop
                        centplusx = fabcenterc + relxcam

                        ; re-init y-loop
                        relycam = (0.0-ir21 + 0.5 -1.0)*fabdyc                               
                        for fj=0,ir[1]-1,incj do begin
                            ; increment relycam
                            relycam = relycam + fabdyc
                            ; update centplusxy to be used in z-loop
                            centplusxy = centplusx + relycam
                            ; reinit z-loop
                            relzcam = (0.0-ir22 + 0.5- 1.0)*fabdzc
                            for fk=0,ir[2]-1,inck do begin                                
                                ; increment relzcam
                                relzcam = relzcam + fabdzc
                                ; update centplusxyz to be used for
                                ; project_cube or a stiple
                                centplusxyz = centplusxy + relzcam
                                if not (((centplusxyz[0]+tempsize[0] gt imglim[0,0]) and (centplusxyz[0]-tempsize[0] lt imglim[1,0])) and $
                                  ((centplusxyz[1]+tempsize[1] gt imglim[0,1]) and (centplusxyz[1]-tempsize[1] lt imglim[1,1])) and $
                                  ((centplusxyz[2]+tempsize[2] gt imglim[0,2]) and (centplusxyz[2]-tempsize[2] lt imglim[1,2]))) then begin
                                    ; cell center is outside camera
                                    centplusxyz = [!values.d_infinity,!values.d_infinity,!values.d_infinity]
                                endif
                                ; no more loops to init
                                
                                ; convert camera space to pixels
                                xind = (centplusxyz[0]-imglim[0,0]) / xcellsize
                                yind = (centplusxyz[1]-imglim[0,1]) / ycellsize
                                zind = (centplusxyz[2]-imglim[0,2]) / zcellsize
                                if sj gt 0.0 then begin
                                                                 
                                    ;if randomu(s) lt (xind-fix(xind)) then xind = fix(xind) else xind = fix(xind+1.0)
                                    ;if randomu(s) lt (yind-fix(yind)) then yind = fix(yind) else yind = fix(yind+1.0)
                                    xind = xind + (randomu(s)-0.5)*sj
                                    yind = yind + (randomu(s)-0.5)*sj
                                    zind = zind + (randomu(s)-0.5)*sj
                                endif
                                xstart = xind - stipsize[0]/2.0
                                xend = xstart + stipsize[0] - 1.0
                                ystart = yind - stipsize[1]/2.0
                                yend = ystart + stipsize[1] - 1.0
                                if usevox then begin
                                    zstart = zind - stipsize[0]/2
                                    zend = zstart + stipsize[0] - 1
                                endif else begin
                                    zstart = 0
                                    zend = 0
                                endelse
                                density = data[fi,fj,fk]
                                ;if density eq 0.0 then begin
                                ;    print,'Density is 0'
                                ;    print,bork
                                ;endif
                                if ((xend gt 0) and (yend gt 0)) and $
                                  ((xstart le ximgsize-1) and (ystart le yimgsize-1)) and density $
                                  then begin
                                ; falls at least somewhat in camera frame
                                ; find if it needs to be clipped
                                    antialias = 0
                                    if antialias eq 0 then begin
                                        dx = xstart - fix(xstart)
                                        dy = ystart - fix(ystart)
                                        xsc = min([round(xstart),0])
                                        xec = max([round(xend)-ximgsize+1,0])
                                        ysc = min([round(ystart),0])
                                        yec = max([round(yend)-yimgsize+1,0])
                                        xstart = round(xstart)
                                        xend = round(xend)
                                        ystart = round(ystart)
                                        yend = round(yend)
                                        if usevox then begin
                                            zsc = min([round(zstart),0])
                                            zec = max([round(zend)-zimgsize+1,0])
                                            if ((-1*xsc lt stipsize[0]) and (stipsize[0]-1 ge xec)) and $
                                              ((-1*ysc lt stipsize[1]) and (stipsize[1]-1 ge yec)) and $
                                              ((-1*zsc lt stipsize[0]) and (stipsize[0]-1 ge zec)) then begin
                                                img[(xstart-xsc):(xend-xec),(ystart-ysc):(yend-yec),(zstart-zsc):(zend-zec)] = $
                                                  img[(xstart-xsc):(xend-xec),(ystart-ysc):(yend-yec),(zstart-zsc):(zend-zec)] $
                                                  + density*Vol1OverVol2*((16.0)^(l))
                                            endif
                                        endif else begin
                                            if ((-1*xsc lt stipsize[0]) and (stipsize[0]-1 ge xec)) and $
                                              ((-1*ysc lt stipsize[1]) and (stipsize[1]-1 ge yec)) then begin
                                                if stipsize[0] ge 2 then begin
                                                    shi = (1.0-dx)*(1.0-dy)*(*stipptrs[l]) + dx*(1.0-dy)*shift((*stipptrs[l]),-1,0)+(1.0-dx)*(dy)*shift((*stipptrs[l]),0,-1)+dx*dy*shift((*stipptrs[l]),-1,-1)
                                                endif else shi = (*stipptrs[l])
                                                img[(xstart-xsc):(xend-xec),(ystart-ysc):(yend-yec)] = $
                                                  img[(xstart-xsc):(xend-xec),(ystart-ysc):(yend-yec)] + $
                                                  shi[(-xsc):(stipsize[0]-xec-1), $
                                                                 (-ysc):(stipsize[1]-yec-1)]*density*VolOverArea
                                            
                                            ;img[(xstart-xsc),(ystart-ysc)] = img[(xstart-xsc),(ystart-ysc)] + density*VolOverArea

                                            endif                                            
                                        endelse
                                    endif else begin ; antialiasing
                                        ;print,'Should not be AA-ing'
                                        ;bork
                                        dx = xind - fix(xind)
                                        dy = yind - fix(yind)
                                        deltx = dx - 0.5
                                        delty = dy - 0.5
                                        ;print,'Dx,Dy = ',string([dx,dy])

                                        ; 0,0 part
                                        xstartp = round(xstart)
                                        xendp = round(xend) 
                                        ystartp = round(ystart)
                                        yendp = round(yend)
                                        xsc =  min([xstartp,0])
                                        xec = max([xendp-ximgsize-1,0])    
                                        ysc = min([ystartp,0])
                                        yec = max([yendp-yimgsize-1,0])
                                        if ((-1*xsc lt stipsize[0]) and (stipsize[0]-1 ge xec)) and $
                                          ((-1*ysc lt stipsize[1]) and (stipsize[1]-1 ge yec)) then begin
                                            img[(xstartp-xsc):(xendp-xec),(ystartp-ysc):(yendp-yec)] = $
                                              img[(xstartp-xsc):(xendp-xec),(ystartp-ysc):(yendp-yec)] + $
                                              (*stipptrs[l])[(-xsc):(stipsize[0]-xec-1), $
                                                             (-ysc):(stipsize[1]-yec-1)]*density*VolOverArea*(1.0-dx)*(1.0-dy)
                                        endif
                                        ;print,'Xi,Yi,Xs,Ys,Xe,Ye',xind,yind,xstartp,ystartp,xendp,yendp
                                        ;print,'Img:  X,X,Y,Y',(xstartp-xsc),(xendp-xec),(ystartp-ysc),(yendp-yec)
                                        ;print,'Stip: X,X,Y,Y',(-xsc),(stipsize[0]-xec-1),(-ysc),(stipsize[1]-yec-1)
                                        
                                        ; -1,0 part
                                        xstartp = round(xstart+deltx)
                                        xendp = round(xend+deltx) 
                                        ystartp = round(ystart)
                                        yendp = round(yend)
                                        xsc =  min([xstartp,0])
                                        xec = max([xendp-ximgsize+1,0])                                       
                                        ysc = min([ystartp,0])
                                        yec = max([yendp-yimgsize+1,0])
                                        if ((-1*xsc lt stipsize[0]) and (stipsize[0]-1 ge xec)) and $
                                          ((-1*ysc lt stipsize[1]) and (stipsize[1]-1 ge yec)) then begin
                                            img[(xstartp-xsc):(xendp-xec),(ystartp-ysc):(yendp-yec)] = $
                                              img[(xstartp-xsc):(xendp-xec),(ystartp-ysc):(yendp-yec)] + $
                                              (*stipptrs[l])[(-xsc):(stipsize[0]-xec-1), $
                                                             (-ysc):(stipsize[1]-yec-1)]*density*VolOverArea*(dx)*(1.0-dy)
                                        endif
                                        ;print,'Xi,Yi,Xs,Ys,Xe,Ye',xind,yind,xstartp,ystartp,xendp,yendp
                                        ;print,'Img:  X,X,Y,Y',(xstartp-xsc),(xendp-xec),(ystartp-ysc),(yendp-yec)
                                        ;print,'Stip: X,X,Y,Y',(-xsc),(stipsize[0]-xec-1),(-ysc),(stipsize[1]-yec-1)

                                        ; 0,-1 part
                                        xstartp = round(xstart)
                                        xendp = round(xend) 
                                        ystartp = round(ystart+delty)
                                        yendp = round(yend+delty)
                                        xsc =  min([xstartp,0])
                                        xec = max([xendp-ximgsize+1,0])                                       
                                        ysc = min([ystartp,0])
                                        yec = max([yendp-yimgsize+1,0])
                                        if ((-1*xsc lt stipsize[0]) and (stipsize[0]-1 ge xec)) and $
                                          ((-1*ysc lt stipsize[1]) and (stipsize[1]-1 ge yec)) then begin
                                            img[(xstartp-xsc):(xendp-xec),(ystartp-ysc):(yendp-yec)] = $
                                              img[(xstartp-xsc):(xendp-xec),(ystartp-ysc):(yendp-yec)] + $
                                              (*stipptrs[l])[(-xsc):(stipsize[0]-xec-1), $
                                                             (-ysc):(stipsize[1]-yec-1)]*density*VolOverArea*(1.0-dx)*(dy)
                                        endif
                                        ;print,'Xi,Yi,Xs,Ys,Xe,Ye',xind,yind,xstartp,ystartp,xendp,yendp
                                        ;print,'Img:  X,X,Y,Y',(xstartp-xsc),(xendp-xec),(ystartp-ysc),(yendp-yec)
                                        ;print,'Stip: X,X,Y,Y',(-xsc),(stipsize[0]-xec-1),(-ysc),(stipsize[1]-yec-1)

                                        ; -1,-1 part
                                        xstartp = round(xstart+deltx)
                                        xendp = round(xend+deltx) 
                                        ystartp = round(ystart+delty)
                                        yendp = round(yend+delty)
                                        xsc =  min([xstartp,0])
                                        xec = max([xendp-ximgsize+1,0])                                       
                                        ysc = min([ystartp,0])
                                        yec = max([yendp-yimgsize+1,0])
                                        if ((-1*xsc lt stipsize[0]) and (stipsize[0]-1 ge xec)) and $
                                          ((-1*ysc lt stipsize[1]) and (stipsize[1]-1 ge yec)) then begin
                                            img[(xstartp-xsc):(xendp-xec),(ystartp-ysc):(yendp-yec)] = $
                                              img[(xstartp-xsc):(xendp-xec),(ystartp-ysc):(yendp-yec)] + $
                                              (*stipptrs[l])[(-xsc):(stipsize[0]-xec-1), $
                                                             (-ysc):(stipsize[1]-yec-1)]*density*VolOverArea*(dx)*(dy)
                                            
                                        endif
                                        ;print,'Xi,Yi,Xs,Ys,Xe,Ye',xind,yind,xstartp,ystartp,xendp,yendp
                                        ;print,'Img:  X,X,Y,Y',(xstartp-xsc),(xendp-xec),(ystartp-ysc),(yendp-yec)
                                        ;print,'Stip: X,X,Y,Y',(-xsc),(stipsize[0]-xec-1),(-ysc),(stipsize[1]-yec-1)

                                    endelse
                                ;barf
                                ;if (xsc ne 0) then begin
                                ;    print,xstart,xend,ystart,yend
                                ;    blargh
                                ;endif
                                    endif ; in frame
                               
                                
                            endfor ;fk                            
                        endfor ;fj
                    endfor ; fi
                    
                    
                    if 0 then begin ; hide project_cube method
                        for fi=ncc.indlo[0], ncc.indhi[0] do begin
                            indfi = fi-ncc.indlo[0]
                            dx = (indfi -  ir20)*tempsize[0]
                            for fj=ncc.indlo[1], ncc.indhi[1] do begin
                                indfj = fj-ncc.indlo[1]
                                dy = (indfj - ir21)*tempsize[1]
                                for fk=ncc.indlo[2], ncc.indhi[2] do begin
                                    indfk = fk-ncc.indlo[2]
                                    dz = (indfk - ir22)*tempsize[2]
                                    celldata = data[indfi,indfj,indfk]
                                    cellcenter = ncc.center $
                                      + [dx,dy,dz]
                                ;print,'Projecting cell with center,density = ',cellcenter,celldata
                                ;project_cube,img,imglim,QSC,celldata,cellcenter,tempsize,subsample
                                endfor
                            endfor
                        endfor
                    endif

                endelse
                
                
            endif else begin
                ;print,'Error in PROJECT_NODE.  Child '+string(chi)+' is of unknown type'
                ;help,ncc
                ;error=1
            endelse
        endif
    endfor

   

endelse


;if  (where(finite(img) ne 1))[0] ne -1 then begin
;    print,'Got NaN!!',bork
;endif

print,'Exiting PROJECT_NODE ',node.level,node.value,(node.xlo+node.xhi)/2.0
if fix(10*randomu(undefinedvar)) eq 6 then begin
    imgm = img
    whi = where(imgm gt 0.0)
    if total(whi) ne -1 then begin
        imm = min(imgm[whi])
    endif else begin
        imm = min(imgm)
    endelse
    whi = where(imgm lt imm)
    if total(whi) ne -1 then begin
        imgm[whi] = imm
    endif
    imgm = alog10(imgm)
    imgm = bytscl(imgm)
    c = intarr(256,3)
    tvlct,c,/get
    imc = [[[imgm]],[[imgm]],[[imgm]]]
    for i=0,2 do imc[*,*,i] = c[imgm[*,*],i]
    tv,imc,true=3
endif

;if  (where(finite(img) ne 1))[0] ne -1 then begin
;    print,'Got NaN!!',bork
;endif

end ; PROJECT_NODE

