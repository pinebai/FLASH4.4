pro slice_node,node,QSC,img,imglim,gridspacing,amrlim

; node is the node to try to take a slice of
; QSC is the projection from sim coords to camera coords
; img is the image to add to
; imglim is the region of the plane represented by img
; gridspacing is the size of a cell on each level
; amrlim is the smallest grid-aligned region of sim space that
; contains imglim


ximgsize = (size(img))[1]
yimgsize = (size(img))[2]
;print,imglim
xpixsize = double((imglim[1,0]-imglim[0,0]))/ double(ximgsize)
ypixsize = double((imglim[1,1]-imglim[0,1]))/ double(yimgsize)

for chin=0,7 do begin
    if ptr_valid(node.children[chin]) then begin
        ncc = *(node.children[chin])
        if n_tags(ncc,/length) eq n_tags({TypeNode},/length) then begin
            if 1 then begin     ; (total( node.xhi lt amrlim[*,0]) + total(node.xlo gt amrlim[*,1])) eq 0 then begin
                slice_node,ncc,QSC,img,imglim,gridspacing,amrlim
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
        if (n_tags(ncc,/length) eq n_tags({TypeFab},/length)) then begin
            ;if 1 lt 2 then print,1
            ; for each fab, perform a cheap test
            ; to see if it is close the the viewing region
            ; if it is, perform a more expensive exact test
            tempsize = double(gridspacing[*,ncc.level])
            ir = (ncc.indhi - ncc.indlo) + [1,1,1]
            fabcenterc = (QSC##[ncc.center,1.0])[0:2]
            maxreach = total((tempsize*ir)/2.0)
            ;print,maxreach
            if 1 then begin
            ;if (abs(fabcenterc[2]) le maxreach) and (total( *(ncc.data) )   ne 0.0) then begin ;((abs(fabcenterc[2]) le maxreach)
                ;print,'Crosses plane'+string(ncc.level)
                if 1 then begin
                ;if ((fabcenterc[1]-maxreach lt imglim[1,1]) and (fabcenterc[1]+maxreach gt imglim[0,1])) then begin
                    if 1 then begin
                    ;if ((fabcenterc[0]-maxreach lt imglim[1,0]) and (fabcenterc[0]+maxreach gt imglim[0,0])) then begin
                        fabdxc = (QSC[0:2,0:2]##[tempsize[0],0,0])
                        fabdyc = (QSC[0:2,0:2]##[0,tempsize[1],0])
                        fabdzc = (QSC[0:2,0:2]##[0,0,tempsize[2]])
                        corners = dblarr(8,3)
                        corners0simspace = ncc.center - [tempsize[0]*ir[0]/2.0, tempsize[1]*ir[1]/2.0 ,tempsize[2]*ir[2]/2.0]   
                        ind = 0
                        above = 0
                        for i=-1,1,2 do begin
                            for j=-1,1,2 do begin
                                for k=-1,1,2 do begin
                                    corners[ind,*] = fabcenterc +i*fabdxc*ir[0]/2.0 +j*fabdyc*ir[1]/2.0 +k*fabdzc*ir[2]/2.0
                                    above = above + (corners[ind,2] gt 0.0)
                                    ind = ind+1
                                endfor ; k loop for corners
                            endfor ; j loop for corners
                        endfor  ; i loop for corners
                        if total(corners[*,2] eq 0.0) eq 4 then above = 0
                        if ((above ne 0) and (above ne 8)) then begin
                            ; fab straddles the plane
                            ;print,([transpose(corners),(dblarr(1,8)+1.0)])
                            QCS = invert(QSC)
                            ;print,'----------------------'
                            ;print,transpose(QCS##transpose([transpose(corners),(dblarr(1,8)+1.0)]))
                            xhilo = 0
                            yhilo = 0
                            for i=0,7 do begin
                                xhilo = xhilo - 1*(corners[i,0] le imglim[0,0]) + (corners[i,0] ge imglim[1,0])
                                yhilo = yhilo - 1*(corners[i,1] le imglim[0,1]) + (corners[i,1] ge imglim[1,1])
                                if (corners[i,0] lt imglim[0,0]) and (corners[i,0] gt imglim[1,0]) then print,bork
                                if (corners[i,1] lt imglim[0,1]) and (corners[i,1] gt imglim[1,1]) then print,bork
                            endfor
                            ;if 1 then begin
                            if (abs(xhilo) ne 8) and (abs(yhilo) ne 8) then begin
                                ; fab interesects viewing window
                                value = ncc.value
                                ximgsize = (size(img))[1]
                                yimgsize = (size(img))[2]
                                xcellsize = double((imglim[1,0]-imglim[0,0]))/ double(ximgsize)
                                ycellsize = double((imglim[1,1]-imglim[0,1]))/ double(yimgsize)
                                ;xmin = min(corners[*,0])
                                ;xmax = max(corners[*,0])
                                ;ymin = min(corners[*,1])
                                ;ymax = max(corners[*,1])
                                
                                if 1 then begin
                                  ; find the points at which the fab
                                  ; edge intersects the plane
                                  ; remake xmax,xmin,ymax,ymin based on
                                  ; this

                                  ; give the order for corner-pairs
                                  ; that make edges
                                    ia = [0,0,0,1,1,2,2,3,4,4,5,6]
                                    ib = [1,2,4,3,5,3,6,7,5,6,7,7]
                                    Xi = 0
                                    Yi = 1
                                    Zi = 2  
                                    ;foo = xmin
                                    ;xmin = xmax
                                    ;xmax = foo
                                    ;foo = ymin
                                    ;ymin = ymax
                                    ;ymax = foo
                                    xmin = !VALUES.F_INFINITY
                                    xmax = -1*xmin
                                    ymin = xmin
                                    ymax = xmax
                                    for e=0,11 do begin
                                        if (corners[ia[e],Zi] gt 0) ne (corners[ib[e],Zi] gt 0)  then begin
                                            p = -1*corners[ib[e],Zi] / ( corners[ia[e],Zi]-corners[ib[e],Zi] )
                                            x = p*corners[ia[e],Xi] + (1-p)*corners[ib[e],Xi]
                                            y = p*corners[ia[e],Yi] + (1-p)*corners[ib[e],Yi]
                                            xmin = min([xmin,x])
                                            xmax = max([xmax,x])
                                            ymin = min([ymin,y])
                                            ymax = max([ymax,y])
                                        endif
                                    endfor                       
                                endif
                                
                                xplo = (xmin - imglim[0,0]) / xcellsize - (1)
                                xphi = (xmax - imglim[0,0]) / xcellsize + (1)
                                yplo = (ymin - imglim[0,1]) / ycellsize - (1)
                                yphi = (ymax - imglim[0,1]) / ycellsize + (1)

                                ;print,'Found FAB:'+strtrim(string(value,xmin,xmax,ymin,ymax),2)
                                ;print,strtrim(string(xplo,xphi,yplo,yphi),2)
                                if (xplo le ximgsize-1) and (xphi ge 0) and (yplo le yimgsize-1) and (yphi ge 0) then begin
                                    ;print,'Level = ',ncc.level
                                    xplo = max([xplo,0])
                                    xphi = min([xphi, ximgsize-1])
                                    yplo = max([yplo,0])
                                    yphi = min([yphi, yimgsize-1])
                                    ;img[xplo:xphi,yplo:yphi] = value
                                    indbase = [xplo*xcellsize + imglim[0,0], yplo*ycellsize + imglim[0,1], 0.0, 1.0]
                                    indbase = (QCS##indbase)[0:2]
                                    indbase = (indbase - corners0simspace)/tempsize
                                    inddx = [xcellsize,       0.0, 0.0]
                                    inddy = [      0.0, ycellsize, 0.0]
                                    inddx = QCS[0:2,0:2]##inddx / tempsize
                                    inddy = QCS[0:2,0:2]##inddy / tempsize
                                    hit = long(0)
                                    miss = long(0)
                                    for xi=xplo,xphi do begin
                                        method = 'by-lines' ; 'by-pixels'  'by-lines' 'by-frames'
                                        if method eq 'by-frames' then begin
                                            if xi gt xplo then xi=xphi+2
                                            xpsz = fix(xphi) - fix(xplo) + 1
                                            ypsz = fix(yphi) - fix(yplo) + 1
                                            inds = dblarr(xpsz,ypsz,3)
                                            for i=0,2 do begin
                                                inds[*,*,i] = indbase[i] + $
                                                  (dindgen(xpsz,ypsz) mod xpsz)*inddx[i] + $
                                                  (dindgen(xpsz,ypsz)  /  xpsz)*inddy[i]
                                                inds[*,*,i] = fix(inds[*,*,i]-0.0)
                                                wh = where( (inds[*,*,i] lt -0.0) or (inds[*,*,i] ge ir[i]-0))
                                                if wh[0] ne -1 then inds[wh,i] = !VALUES.F_NAN
                                            endfor
                                            lininds = inds[*,*,0] + inds[*,*,1]*ir[0] + inds[*,*,2]*ir[0]*ir[1]
                                            wh = where( ((*(ncc.data))[lininds]) eq 0.0)
                                            if wh[0] ne -1 then lininds[wh] = !VALUES.F_NAN
                                            goodmask = finite(lininds)
                                            good = where(goodmask)
                                            if good[0] ne -1 then begin
                                                help,img
                                                print,xplo,xphi,yplo,yphi
                                                imgsub = img[xplo:xphi,yplo:yphi]
                                                imgsub[good] = ((*(ncc.data))[lininds]*goodmask)[good]
                                                img[xplo:xphi,yplo:yphi] = imgsub
                                                ;(img[xplo:xphi,yplo:xphi])[good] = ((*(ncc.data))[lininds]*goodmask)[good]
                                            endif
                                        endif
                                        if method eq 'by-lines' then begin
                                            ypsz = fix(yphi) - fix(yplo) + 1
                                            inds = dblarr(ypsz,3)
                                            for i=0,2 do begin
                                                inds[*,i] = indbase[i] + ((xi-xplo)*inddx)[i] + ((indgen(ypsz)+0*yplo)*inddy[i])
                                                inds[*,i] = fix(inds[*,i] - 0.0)
                                                wh = where( (inds[*,i] lt -0.0) or (inds[*,i] ge ir[i]-0))
                                                if wh[0] ne -1 then inds[wh,i] = !VALUES.F_NAN
                                            endfor
                                            ;if ncc.level gt 1 then print,bork
                                            lininds = inds[*,0] + inds[*,1]*ir[0] + inds[*,2]*ir[0]*ir[1]
                                            ;print,bork
                                            wh = where( ((*(ncc.data))[lininds]) eq 0.0)
                                            if wh[0] ne -1 then lininds[wh] = !VALUES.F_NAN
                                            goodmask = finite(lininds)
                                            good = where(goodmask)
                                            ;goodmask = goodmask and ((*(ncc.data))[lininds]*goodmask)
                                            ;good = where(goodmask)
                                            if good[0] ne -1 then img[xi,good+yplo] = ((*(ncc.data))[lininds]*goodmask)[good]
                                            ;if good[0] ne -1 then begin
                                            ;    plot,good+yplo
                                            ;endif else begin
                                            ;endelse
                                        endif
                                        ;if good[0] ne -1 then print,bork
                                        ;This is correct except that if inds[*,0] or inds [*,1] or inds[*,2] exceed their respective bounds, then lininds must make sure to not reverence it,  try linds[offending value] = -1
                                        ;print,bork
                                        ;lininds = inds[*,0] + inds[*,1]*ir[0] + inds[*,2]*ir[0]*ir[1]         
                                        ;wherebad = where( (inds lt 0.0 ) or (inds ge ir[0]*ir[1]*ir[2]-1))
                                        ;vals = (*(ncc.data))[lininds]
                                        ;img[xi,yplo:yphi] = vals

                                        ;inds = linearized indices of data that correspond to the pixels
                                        ;some inds are out of range
                                        ;goodmask = (inds ge minimum_legal_value) and (inds le maximum_legal_value)
                                        ;goodmask is a mask of which inds are valid
                                        ;data[inds]*goodmask is values to assign to the pixels, or 0 if invalid
                                        ;where(goodmask) is the addrs of the image that match parts of the FAB
                                        ;where(goodmask) is also the arrds of inds that are valid
                                        ;dest[where(goodmask)] = (data[inds]*goodmask)[where(goodmask)]
                                        ; ABOVE sets only pixels matching the fab to the parts of the fab they match
                                        ;print,bork
                                        if method eq 'by-pixels' then begin
                                        for yi=yplo,yphi do begin
                                            ; convert pixel location to simspace.
                                ; if in fab, set pixel value to cell
; value
                                            ; old method
                                            ;planecoord = [xi*xcellsize + imglim[0,0], yi*ycellsize + imglim[0,1], 0.0, 1.0]
                                            ;simcoord = (QCS##planecoord)[0:2]                                
                                ;index = (simcoorrd - corners0simspace)/tempsize
                                            index = indbase + (xi-xplo)*inddx + (yi-yplo)*inddy
                                            ;print,bork

                                            ;print,corners0simspace,index
                                            ; index = (simcoord - corners[0,*])/tempsize

                                            if (total(index gt -1.0) eq 3) and (total(index lt ir) eq 3) then begin   
                                                value = (*(ncc.data))[index[0],index[1],index[2]]
                                                if value then begin
                                                    img[xi,yi] = value
                                                    hit = hit + long(1)
                                                endif else begin
                                                    miss = miss + long(1)
                                                endelse
                                            endif else begin  
                                                miss = miss + long(1)
                                            endelse
                                        endfor ; y-loop                              
                                    endif ; if by-pixels

                                endfor ; x-loop
                                    ;print,'Hit, miss, pct',hit,miss,double(hit)/double(hit+miss)
                                    ;if hit eq long(0) then print,bork
                                endif
                                ; foreach pixel in
                                ; img[xmin:xmax,ymin:ymax] test if it in in the cube
                                
                            endif
                        endif
                    endif       ; x-test for being close to window
                endif           ; y-test for being close to window
            endif               ; fab is close to plane
        endif                   ;is FAB
    endif                       ;if valid pointer
endfor
end

