pro node_describe,node
  print,'Node:'
  print,node
  for i=0 , 7 do begin
      if ptr_valid(node.children[i]) then begin
          if n_tags( *(node.children[i])) eq n_tags({TypeNode}) then begin
              node_describe,*(node.children[i])
          endif else if n_tags( *(node.children[i])) eq n_tags({TypeFab}) then begin
              fab_describe,*(node.children[i])
          endif
      endif
  endfor
end
