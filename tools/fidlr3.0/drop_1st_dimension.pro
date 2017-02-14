;--------------------------------------------------------------
; drop_1st_dimension.pro
;
; Take an array 'arr' of n dimensions whose first dimension is
; assumed to be degenerate, that is, of size=1, and return an
; array of n-1 dimensions where the degenerate first dimension
; has been removed.
;
; drop_1st_dimension offers improved control over a direct use
; of idl's built-in "reform" function, which always removes all
; dimensions of size=1, a behavior which for flash will cause
; errors in simulations done on only one block.
; 
; The function raises an error if the first dimension of 'arr'
; is not degenerate.
;--------------------------------------------------------------
function drop_1st_dimension, arr
  return, reform(arr, (size(arr, /DIMENSIONS))[1:*])
end
