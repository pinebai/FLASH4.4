function amr_line,amr=amr,tree=tree,vector=vector,point=point,range=range,maxlevel=maxlevel,numpoints=numpoints,_extra=_extra

; Copyright Daniel Gies, 2003
; License:  Public Domain
; No warranty, either express or implied.  Use of this software means
; you assume all liability arising from its use.


; Program to extract data which falls along a line through an AMR or
; amr oct-tree data structure

; Input Variables
; amr : AMR dataset to be examined
; tree : AMR oct-tree to be examined
; vector:  [x,y,z] vector along which to extract data
; point: [x,y,z] point the vector should originate at
; range: [-a,+b] range along vector from point to examine
; maxlevel: sets max level to extract data from
; numpoints: length of vector to be returned
; Either amr or tree must be specified.  If both are given, amr is
; used

; Output variable is a dblarr() with length determined by maxlevel and
; range, or by numpoints if specified

MNP = 65536
; Hard limit on maximum number of points
ndim=3

if keyword_set(amr) then begin
    ; work on amr object
    if not keyword_set(vector) then begin
        vector = dblarr(ndim)        
        vector[0] = 1.0
    endif
    ; normalize vector
    vector = vector / sqrt(total(vector*vector))
    if not keyword_set(point) then begin
        point = (amr.boxmin+amr.boxmax)/2.0
    endif
    if not keyword_set(range) then begin
        range=[amr.boxmin[0],amr.boxmax[0]]
        ; find by intersection of vector with amr domain
        ; incomplete
    endif
    if not keyword_set(maxlevel) then begin
        maxlevel = amr.maxlevel
    endif
    if maxlevel gt amr.maxlevel then maxlevel=amr.maxlevel
    if not keyword_set(numpoints) then begin
        numpoints = double(range[1]-range[0])/amr.gridspacing[0,maxlevel]
        numpoints = long(numpoints)
    endif
    ; clip numpoints to a sensible value
    if numpoints gt MNP then numpoints = MNP
    if numpoints lt 2 then numpoints = long(2)
    ; create result vector
    result = dblarr(numpoints)
    
    ; start and end points of the search line
    p1 = point + vector*range[0]
    p2 = point + vector*range[1]
    ; amount to step through for each point in the result
    dp = double(range[1]-range[0])/double(numpoints-1)
    if 1 then begin
        for i=long(0),numpoints-1 do begin
            xyz = p1+i*dp
            coord = coord_to_fab(xyz,amr)
            result[i] = coord_to_val(xyz,amr)
            ;if coord[0] ne -1 then begin
            ;    ml=0
            ;    if (size(coord))[0] eq 2 then begin
            ;        ml = (size(coord))[2]-1
            ;        coord = coord[*,ml]
            ;        ;help,ml
            ;    endif
            ;    l = coord[0]
            ;    f = coord[1]
            ;    ix = coord[2]
            ;    iy = coord[3]
            ;    iz = coord[4]
            ;    result[i] = (*((*amr.levels[l].fabptr)[f].dataptr))[ix,iy,iz]
            ;endif else begin
            ;    result[i] = !values.f_nan
            ;    ;print,xyz
            ;endelse
        endfor
    endif
    if 0 then begin
    for i=long(0),numpoints-1 do begin
        inc = 0
        pointi = p1+i*dp
        ; find finest grid that contains pointi
        fabid=[-1,-1]
        for l=maxlevel,0,-1 do begin
            for f=0,amr.levels[l].nfab-1 do begin
                if total(pointi ge (*amr.levels[l].fabptr)[f].xlo) eq ndim $
                  and total(pointi le (*amr.levels[l].fabptr)[f].xhi) eq ndim then begin
                    ; pointi is in fab f on level l
                    ; record the fabid
                    fabid = [l,f]
                    ; cause inner 2 loops to terminate
                    f = amr.levels[l].nfab
                    l = -1
                endif
            endfor
        endfor
        ; if grid does not exist, result[i] = NaN
        ; else iterate on that grid until hit boundary
        if fabid[1] eq -1 then begin
            result[i] = !values.f_nan
        endif else begin
            fp = (*amr.levels[l].fabptr)[f]
            xyzmin = ((*fp).idxlo-((*fp).idxhi-(*fp).idxlo)/2.0) * amr.gridspacing[*,l]
            xyzmax = ((*fp).idxhi-((*fp).idxhi-(*fp).idxlo)/2.0) * amr.gridspacing[*,l]
            print,'XYZminmax:',xyzmin,xyzmax
            for j=pointi,numpoints-1 do begin
                                ; if xi,yi,zi out of range, j=numpoints
                                ; else result[j] = data[xi,yi,zi]
                xyz = p1+j*dp
                if (total(xyz gt xyzmax) ne 0) or (total(xyz lt xyzmin) ne 0) then begin
                    inc = j-pointi
                    j = numpoints
                endif else begin
                    
                endelse
            endfor
            ; increase i by number of iterations in inner loop
            i = i+inc
        endelse
    endfor
endif


endif else if keyword_set(tree) then begin
    ; work on amr oct-tree
endif
return,result
end
