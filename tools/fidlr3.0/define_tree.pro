;===============================================================================
; define a data structure for the tree in double or single precision
; Arguments:  tree (out)  -- the tree structure
;             precision (in) -- 1 return the variables in single precision
;			     -- 2 return the variables in double precision
;
;;==============================================================================

pro define_tree, precision, nfaces, nchild, ndim, TREE=tree

; if no precision is passed, presume double
IF (n_elements(precision) EQ 0) THEN BEGIN
    precision = 2
    print,"No precision requested, returning double precision by default"
ENDIF 
IF (n_elements(nfaces) EQ 0) OR (n_elements(nchild) EQ 0) THEN BEGIN
    print," Failure in define_tree, must call with nfaces and nchild"
ENDIF 
IF (n_elements(ndim) EQ 0) THEN BEGIN
    print, "Number of dimensions not requested, assuming 2"
    ndim = 2
ENDIF 

;------------------------------------------------------------------------------
;   setup the structures to pass the data
;----------------------------------------------------------------------------

IF ( precision EQ 1) then begin
    tree = {lrefine:0l, $
            nodeType:0l, $
            processorNumber:0l, $
            gid:lonarr(nfaces+1+nchild), $
            coord:fltarr(ndim), $
            size:fltarr(ndim), $
            bndBox:fltarr(2,ndim)}
ENDIF   
IF ( precision EQ 2) THEN begin
    tree = {lrefine:0l, $
            nodeType:0l, $
            processorNumber:0l, $
            gid:lonarr(nfaces+1+nchild), $
            coord:dblarr(ndim), $
            size:dblarr(ndim), $
            bndBox:dblarr(2,ndim)}
ENDIF 

; following line needs to be called in main program
;tree = replicate(tree,totBlocks)

END   ; of define_tree.pro    

