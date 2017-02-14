pro node_addfab,node,fab,gridspacing,quiet=quiet

;anode = {TypeNode,children:ptrarr(8), xlo:dblarr(3), xhi:dblarr(3),
;level:fix(0), fabsxlo:dblarr(3), fabsxhi:dblarr(3)}
;acell = {TypeCell,value:double(0.0)}





; node:Node is the node to have fab added to
; fab:Fab os the Fab to be glommed onto the tree
; center:fltarr(ndim) is the center of the fab
; cellsize:fltarr(ndim) is the size of a cell in the fab
; level:int is the level the fab lives on
error=0
center = fab.center
if total(center lt node.xlo) then begin
    print,'Error! Fab falls outside node! (low)'
    print,'Node: '+string(node.xlo)+' - '+string(node.xhi)
    print,'Cell: '+string(center)
    error=1
endif 
if total(center gt node.xhi) then begin
    print,'Error! Fab falls outside node! (high)' 
    print,'Node: '+string(node.xlo)+' - '+string(node.xhi)
    print,'Cell: '+string(center)
    error=1
endif
if error ne 0 then begin
    return
endif

nodec = (node.xlo+node.xhi)/2.0
; find out which octant it goes into
octant = 0
;center = fab.center
if center[0] lt nodec[0] then begin
    octant = octant+0
endif else begin 
    octant = octant+1
endelse

if center[1] lt nodec[1] then begin
    octant = octant+0
endif else begin
    octant = octant+2
endelse

if center[2] lt nodec[2] then begin
    octant = octant+0
endif else begin 
    octant = octant+4
endelse
if not keyword_set(quiet) then begin
    print,'Node at depth '+string(node.level)+', octant='+string(octant)+' level='+string(node.level)
endif
;if ptr_valid(node.children[octant]) then begin
;    print,'Ntags='+string(n_tags( *(node.children[octant])))
;endif else begin
;    print,'Was empty, is nod filled:'
;    fab_describe,fab
;endelse

; if octant empty, make pointer to Fab type
; if octant full, and points to node, recurse function on child node
; if octant full and points to fab of same level , make node in child and fill with
; 2 fabs
; but if full with fab of lower level, kill fab and replace with new
; node pointing to new fab
if ptr_valid(node.children[octant]) then begin
    ; test if pointer to node or fab
    ; WARNING: Will crash if TypeFab and
    ; TypeNode take the same # of bytes.  I don't know how to test for
    ; structure type in IDL
    if n_tags( *(node.children[octant]),/length) eq n_tags({TypeFab},/length) then begin
        ; pointer to a fab
        ; replace fab with a new node and
        ; call node_addfab on it with both fabs
        oldfab = *(node.children[octant])
        if oldfab.level lt fab.level then begin            
                ; test for same volume, center, and aspect ratios
            if 0 then begin; 2 fabs on differenet levels cover each other exactly
                                ;if (total(fab.center eq oldfab.center) eq 3) and $
                                ;  (fab.dv eq oldfab.dv ) and $
                                ;  ( stddev((fab.indhi-fab.indlo) / (oldfab.indhi-oldfab.indlo)) eq 0 ) $
                                ; if new fab completely covers oldfab,
                                ; replace oldfab
                if not keyword_set(quiet) then begin
                    print,'Replacing octant '+string(octant)+' level '+string(oldfab.level)+' fab with level '+string(fab.level)+' fab.'
                endif
                node.value = node.value - oldfab.value
                node.value = node.value + fab.value
                ptrtmp = node.children[octant]
                ptrd = (*(node.children[octant])).data
                ptr_free,ptrd
                ptr_free,ptrtmp
                node.children[octant] = ptr_new(fab,/no_copy)
            endif else begin
                                ; else replace with node that points
                                ; to both fabs, and zero out parts of
                                ; oldfab covered by newfab
                rr = gridspacing[*,fab.level] / gridspacing[*,oldfab.level]
                olddatp = oldfab.data
                origoldfabvalue = oldfab.value
                                ;for i=fab.indlo[0],fab.indhi[0] do begin
                                ;    for j=fab.indlo[1],fab.indhi[1] do begin
                                ;        for k=fab.indlo[2],fab.indhi[2] do begin
                                ;            (*olddatp)[i*rr,j*rr,k*rr] = 0.0
                                ;        endfor
                                ;    endfor
                                ;endfor ; i
                                ; recompute oldfab.value 
                oldfab.value = total(*(oldfab.data)) * oldfab.dv
                node.value = node.value - origoldfabvalue + oldfab.value + fab.value
                ptr_free,node.children[octant]
                nxlo = node.xlo
                nxhi = node.xhi        
                if center[0] lt nodec[0] then begin
                    nxhi[0] = nodec[0]
                endif else begin 
                    nxlo[0] = nodec[0]
                endelse                    
                if center[1] lt nodec[1] then begin
                    nxhi[1] = nodec[1]
                endif else begin
                    nxlo[1] = nodec[1]
                endelse
                if center[2] lt nodec[2] then begin
                    nxhi[2] = nodec[2]
                endif else begin 
                    nxlo[2] = nodec[2]
                endelse
                node.children[octant] = $
                  ptr_new({TypeNode,children:ptrarr(8), $
                           xlo:nxlo, xhi:nxhi, level:fix(node.level + 1), value:double(0.0), fabsxlo:dblarr(3), fabsxhi:dblarr(3)})                    
                if oldfab.level ne fab.level then begin
                    if not keyword_set(quiet) then begin
                        print,'Warning, making Node from 2 Fabs with different levels, ',fab.level,oldfab.level
                    endif
                endif
                if not keyword_set(quiet) then begin
                    print,'Replacing octant '+string(octant)+' level '+string(oldfab.level)+' fab with node and level '+string(fab.level)+' fab'
                endif
                nco = *(node.children[octant])
                if not keyword_set(quiet) then begin
                    print,nco
                    print,oldfab
                    print,fab
                endif
                node_addfab,nco,oldfab,gridspacing,quiet=quiet
                *(node.children[octant]) = nco
                
                nco = *(node.children[octant])
                node_addfab,nco,   fab,gridspacing,quiet=quiet
                *(node.children[octant]) = nco
            endelse
        endif else begin
                                ; if oldfab is same level, replace
                                ; oldfab with node that points to fab and oldfab
            oldfab = *(node.children[octant])
            if n_tags(oldfab,/length) ne n_tags({TypeFab},/length) then begin
                print,'Error! oldfab not a fab!'
                blargh
            endif
            ptr_free,node.children[octant]
            ;node.children[octant] = $
            ;  ptr_new({TypeNode,children:ptrarr(8), $
            ;           xlo:(node.xlo+center)/2.0, $
            ;           xhi:(node.xhi+center)/2.0,$
            ;           level:fix(node.level + 1)})
            nxlo = node.xlo
            nxhi = node.xhi
            ;;;;
            if center[0] lt nodec[0] then begin
                nxhi[0] = nodec[0]
            endif else begin 
                nxlo[0] = nodec[0]
            endelse
            
            if center[1] lt nodec[1] then begin
                nxhi[1] = nodec[1]
            endif else begin
                nxlo[1] = nodec[1]
            endelse
            
            if center[2] lt nodec[2] then begin
                nxhi[2] = nodec[2]
            endif else begin 
                nxlo[2] = nodec[2]
            endelse
            ;;;;
            node.children[octant] = $
              ptr_new({TypeNode,children:ptrarr(8), $
                       xlo:nxlo, xhi:nxhi, level:fix(node.level + 1), value:double(0.0), fabsxlo:dblarr(3), fabsxhi:dblarr(3) })
            
            if oldfab.level ne fab.level then begin
                if not keyword_set(quiet) then begin
                    print,'Warning, making Node from 2 Fabs with different levels, ',fab.level,oldfab.level
                endif
            endif
            if not keyword_set(quiet) then begin
                print,'Replacing octant '+string(octant)+' level '+string(oldfab.level)+' fab with node and level '+string(fab.level)+' fab'
            endif
            nco = *(node.children[octant])
            if not keyword_set(quiet) then begin
                print,nco
                print,oldfab
                print,fab
            endif
                                ;blargh
            node_addfab,nco,oldfab,gridspacing,quiet=quiet
            *(node.children[octant]) = nco
            
            node.value = node.value + fab.value*fab.dv
            nco = *(node.children[octant])
            node_addfab,nco,   fab,gridspacing,quiet=quiet
            *(node.children[octant]) = nco
        endelse
    endif else begin
                                ; pointer to a node
                                ; call node_addfab on it
        if not keyword_set(quiet) then begin
            print,'Recursing on octant '+string(octant)
        endif
        nco = *(node.children[octant])
        childvaluebefore = nco.value      
        node_addfab,nco,fab,gridspacing,quiet=quiet
        node.value = node.value - childvaluebefore + nco.value
        *(node.children[octant]) = nco
    endelse                     ; does not point to a fab
endif else begin                ; octant is empty
    if not keyword_set(quiet) then begin
        print,'Making empty octant'+string(octant)+' point to level '+string(fab.level)+' fab.'
    endif
    node.value = node.value + fab.value*fab.dv
    node.children[octant] = ptr_new(fab,/no_copy)
                                ;print,'is valid? ',ptr_valid(node.children[octant])
endelse                         ; children[octant] is null pointer

;print,node
if total(ptr_valid(node.children)) eq 0 then begin
print,'Error! No children!  Borking!'
foo,bar
endif
end                             ;node_addfab
