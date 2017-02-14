pro compare, FILE1 = file1, FILE2 = file2

; crude little file comparison tool
; NOTE!  Morton ordering is different in flash2 versus flash3
;  So this ditty is likely to fail.  
;  Use diff.pro instead


if (n_elements(file1) EQ 0) then begin
    print, 'ERROR, must specify FILE1'
    return
endif

if (n_elements(file2) EQ 0) then begin
    print, 'ERROR, must specify FILE2'
    return
endif

read_amr, file1, PARAMETERS=params1, TREE=tree1, DATA=data1, $
  STORED_VARS = vars1
read_amr, file2, PARAMETERS=params2, TREE=tree2, DATA=data2, $
  STORED_VARS = vars2

; check to make sure that the number of variables agree
if (params1.nvar NE params2.nvar) then begin
    print, 'ERROR: number of variables do not agree'
    return
endif

; check to see if the block structure agrees -- eventually this should
; tell you where they disagree and print out the different numbers
if (params1.totBlocks NE params1.totBlocks) then begin
    print, 'ERROR: number of blocks does not agree'
    return
endif

refine1 = tree1[*].lrefine
refine2 = tree2[*].lrefine
refinediff = refine1-refine2
diff = max(abs(refinediff))
if (diff NE 0) then begin
    print, 'ERROR: refinement levels do not agree'
    print, 'Maximum refinement for file 1 is ',max(refine1)
    print, 'Maximum refinement for file 2 is ',max(refine2)
    print, 'Refinement structure for file 1 is ',refine1
    print, 'Refinement structure for file 2 is ',refine2
    print, 'Difference structure is ',refinediff
    ; refinement may not be a fatal flaw... continue
    ; return
endif

node1 = tree1[*].nodetype
node2 = tree2[*].nodetype
diff = max(abs(node1 - node2))
if (diff NE 0) then begin
    print, 'ERROR: nodetype does not agree'
    print, 'Nodes of file 1 are ',node1
    print, 'Nodes of file 2 are ',node2
    ; again, probably not a fatal flaw....
    ;return
endif

coord1 = tree1[*].coord
coord2 = tree2[*].coord
diff = max(abs(coord1-coord2))
if (diff) then begin
    print, 'ERROR: coords do not agree'
    print, 'Coordinates of file 1 are ',coord1
    print, 'Coordinates of file 2 are ',coord2
    ;return
endif

diff = max(abs(tree1[*].size - tree2[*].size))
if (diff) then begin
    print, 'ERROR: size does not agree'
    return
endif


; check the individual variables
for i = 0, params1.nvar-1 do begin
    if (vars1[i] NE vars2[i]) then begin
        print, 'ERROR: variable names do not match'
        print, 'file1: ', vars1[i], ' file2: ', vars2[i]
        return
    endif

    diff = max(abs(data1[i,*,*,*,*] - data2[i,*,*,*,*]))
    if (diff NE 0.d0) then begin
        print, 'ERROR: differences in var: ', vars1[i]
        var1 = data1[i,*,*,*,*]
        var2 = data2[i,*,*,*,*]
        i_nonzero = where(var1 NE 0)

        print, '       relative error = ', $
          max(abs((var1[i_nonzero] - var2[i_nonzero])/var1[i_nonzero]))
        index = where(max(abs((var1[i_nonzero] - var2[i_nonzero])/var1[i_nonzero])) EQ $
                      abs((var1[i_nonzero] - var2[i_nonzero])/var1[i_nonzero]))
        print, 'loc = ', var1[i_nonzero[index[0]]], var2[i_nonzero[index[0]]]
        print, '       absolute error = ', $
          max(abs((var1 - var2)))
        index = where(max(abs(var1 - var2)) EQ abs(var1 - var2))
        print, 'loc = ', var1[index[0]], var2[index[0]]
    endif

endfor

end


