
function amr_tree,amr,maxlevel=maxlevel,free=free,quiet=quiet
; Takes an AMR data structure and rearranges all cells into
; an octtree using only the finest data available at each level
; if maxlevel is set, do not use any data finer than maxlevel
; if free is true, free data in the amr structure as you go along

anode = {TypeNode, $
         children:ptrarr(8), $
         xlo:dblarr(3), $
         xhi:dblarr(3), $
         level:fix(0), $
         value:double(0.0), $
         fabsxlo:dblarr(3), $
         fabsxhi:dblarr(3) $
        }
acell = {TypeCell,value:double(0.0)}
afab = {TypeFab, $
        center:dblarr(3), $
        indlo:lonarr(3), $
        indhi:lonarr(3),  $
        level:fix(0),$
        value:double(0.0),$
        dv:double(0.0), $ 
        data:ptr_new() $
       }


if amr.ndim ne 3 then begin
  print,'Error: amr_tree only supports 3-D data at this time'
  return,0
endif

name = amr.name
componentName = amr.componentName
version = amr.version
ndim = amr.ndim
time = amr.time
if not keyword_set(maxlevel) then begin 
    maxlevel = amr.maxlevel
endif else begin 
    if maxlevel gt amr.maxlevel then maxlevel = amr.maxlevel
endelse
boxmin = amr.boxmin
boxmax = amr.boxmax
refratio = amr.refratio
idxlo = amr.idxlo
idxhi = amr.idxhi
idxtype = amr.idxtype
periodicity = amr.periodicity
levelsteps = amr.levelsteps
gridspacing = amr.gridspacing
coordtype = amr.coordtype
node = ptr_new(Node)

tree = {TypeAMRTree,name:amr.name, componentName:amr.componentName, $
           version:amr.version, ndim:amr.ndim, $
           time:amr.time, maxlevel:maxlevel, $
           boxmin:amr.boxmin, boxmax:amr.boxmax, $
           refratio:amr.refratio, $ ; contains extra info if (maxlevel ne amr.maxlevel)
           gridspacing:amr.gridspacing, $
           maxindex:max(amr.idxhi[*,amr.maxlevel-1]), $
           srange:dblarr(2), $
           rootnode:{TypeNode}}

tree.rootnode = {TypeNode}
tree.rootnode.xlo = tree.boxmin
tree.rootnode.xhi = tree.boxmax
tree.rootnode.level = fix(0)
tree.rootnode.value = double(0.0)
scalarmin = min_amr(amr)
scalarmax = max_amr(amr)
tree.srange =[scalarmin,scalarmax]

; zero-out redundant data
; for each level higher that 0, check the next lowest level for
; overlap
zerooutoverlap = 1
if zerooutoverlap ne 0 then begin
    for levi=1, maxlevel do begin
        for fabzi=0, amr.levels[levi].nfab-1 do begin
                                ; check if fab imtersects any lower-level fabs
                                ; if it does, zero-out those parts
                                ; convert indlo and indhi into ranges in next-lower level
            indil = ((*(amr.levels[levi].fabptr))[fabzi]).idxlo / amr.refratio[levi-1]
            indih = ((*(amr.levels[levi].fabptr))[fabzi]).idxhi / amr.refratio[levi-1]
            levj = levi-1
            
            for fabzj=0, amr.levels[levj].nfab-1 do begin
                                ; (*(amr.levels[levi].fabptr))[fabzi] 
                                ; (*(amr.levels[levj].fabptr))[fabzj] 
                indjl = ((*(amr.levels[levj].fabptr))[fabzj]).idxlo
                indjh = ((*(amr.levels[levj].fabptr))[fabzj]).idxhi
                if (total(indil le indjh) eq 3) and (total(indih ge indjl) eq 3) then begin
                    if not keyword_set(quiet) then begin
                        print,'Zeroing fab at level'+string(levj)
                    endif
                                ; fabs intersect
                                ; loop on indices of intersection and
                                ; zero-out lower level fab where it intersects
                    indstart = lonarr(3)
                    indend = lonarr(3)
                    indstart[0] = max([indil[0], indjl[0]])
                    indstart[1] = max([indil[1], indjl[1]])
                    indstart[2] = max([indil[2], indjl[2]])
                    indend[0] = min([indih[0], indjh[0]])
                    indend[1] = min([indih[1], indjh[1]])
                    indend[2] = min([indih[2], indjh[2]])
                    indstart = indstart - indjl
                    indend = indend - indjl
                    if 1 then begin
                        (*(((*(amr.levels[levj].fabptr))[fabzj]).dataptr))$
                          [indstart[0]:indend[0],$
                           indstart[1]:indend[1],$
                           indstart[2]:indend[2]]$
                          = 0.0
                    endif
                endif
            endfor
        endfor
    endfor
endif

; for each level <= maxlevel
;   for each fab
;     node_addfab,rootnode,fab
;     if free, free amr memory

for lev=0, maxlevel do begin
  for fabi=0, amr.levels[lev].nfab-1 do begin
      fab = {TypeFab}
      amrfp = *(amr.levels[lev].fabptr)
      fab.indlo = amrfp[fabi].idxlo
      fab.indhi = amrfp[fabi].idxhi
      fab.center = (amrfp[fabi].xlo + amrfp[fabi].xhi)/2.0
      fab.level = lev
      fab.dv = double(tree.gridspacing[0,lev])*double(tree.gridspacing[1,lev])*double(tree.gridspacing[2,lev])
      fab.value = total(*(amrfp[fabi].dataptr)) * fab.dv
      if not keyword_set(quiet) then begin
          print,'Value= '+string(fab.value)
      endif
      if free then begin         
          fab.data = amrfp[fabi].dataptr
          amrfp[fabi].dataptr = ptr_new() ; make a null pointer
      endif else begin
          fab.data = ptr_new( *(amrfp[fabi].dataptr)) ; make a copy
      endelse
      if not keyword_set(quiet) then begin
          print,'Calling node_addfab on a level '+string(lev)+' fab'
          print,tree.rootnode
          print,fab
      endif
      rootnode = tree.rootnode
      node_addfab,rootnode,fab,tree.gridspacing,quiet=quiet
      tree.rootnode = rootnode
  endfor ; fab
endfor ; level
return,tree
end ; amr_tree



;function cell_describe,cell
;  print,'Value = '+cell.value+', Level = '+cell.level
;end




;function tree_addfab,tree,fab
;error=0
;if size(where(cellcenter lt node.xlo),/n_elements) then begin
;    print,'Error! Fab falls outside node! (low)'
;    error=1
;endif 
;if size(where(cellcenter gt node.xhi),/n_elements) then begin
;    print,'Error! Fab falls outside node! (high)'
;    error=1
;endif
;if error then begin
;    print,'Error in tree_addfab !  Exiting'
;    node_describe,tree.rootnode
;    fab_describe,fab
;endif else begin
;    node_addfab,tree.rootnode,fab
;endelse
;end ; tree_addfab
