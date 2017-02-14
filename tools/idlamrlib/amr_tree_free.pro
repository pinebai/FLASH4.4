pro amr_tree_free,node

for i=0,7 do begin
  if ptr_valid(node.children[i]) then begin
      if n_tags( *(node.children[octant])) ne 4 then begin ; a fab
          if ptr_valid( *(node.children[octant]).data ) then begin
              ptr_free(*(node.children[octant]).data)
          endif
      endif else begin ; a node
          amr_tree_free,*(node.children[octant])
      endelse
      ptr_free(node.children[i])
  endif
endfor
