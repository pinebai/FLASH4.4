pro diff3, filename1,varname1,filename2,varname2,debug,OUTPUT=strOutputType, NOINTERP=nointerpolation

; some variables will be saved
common save_ptr, tree_ptr, data_ptr, params_ptr, xmerge_ptr, ymerge_ptr, zmerge_ptr
common save, current_block, current_var
common variables, varnames, native_num

if N_ELEMENTS(strOutputType) EQ 0 then begin
    outputType = 0
endif else begin
    if (strOutputType EQ "ps") then begin 
        outputType=1
    endif else begin
        outputType = 0
    end
endelse

nointerp = Keyword_Set(nointerpolation)

; get the xflash directory from the xflash_dir environmental variable
xflash_dir = get_xflash_path()
; find the current path
spawn, 'pwd', current_path
current_path = current_path[0]
print, '... current path is ', current_path
print, '... XFLASH_DIR is ', xflash_dir


; read in the file, and make sure it is a real FLASH file
itype = determine_file_type(filename1)
itype2 = determine_file_type(filename2)

if (itype EQ -1 OR itype2 EQ -1) then begin
    print, 'ERROR: invalid file'
    return
endif


ndim = determine_file_dimensionality(filename1)
ndim2 = determine_file_dimensionality(filename2)

if (ndim NE ndim2) then begin
    print, 'ERROR: files have different dimensionality'
    return
endif

tree = {lrefine:0l, $
        nodeType:0l, $
        gid:lonarr(2*ndim+1+2^ndim), $
        coord:fltarr(ndim), $
        size:fltarr(ndim), $
        bndBox:fltarr(2,ndim)}

tree_ptr1 = ptr_new(tree)
tree_ptr2 = ptr_new(tree)

params = {totBlocks:0, $
          corners:0, $
          ndim:0, $
          nvar:0, $
          nxb:0, $
          nyb:0, $
          nzb:0, $
          ntopx:1, $
          ntopy:1, $
          ntopz:1, $
          time:0.0, $
          dt:0.0}
params_ptr1 = ptr_new(params)
params_ptr2 = ptr_new(params)
data_ptr1 = ptr_new(0.0)
data_ptr2 = ptr_new(0.0)

xflash_defaults, 'Generic', $
  CTR_LIM=ctr_lim, ISWP=iswp, $
  CONTOURS=contours_def, VECTOR=vector_def, $
  NBINS = nbins, HIST_SCALE = max_scale
                        
; save the contours
contours = contours_def
vector = vector_def
                
; we need to save ctr_lim, so other events can access it
pblm_ctr_lim = ctr_lim
orientation = iswp

current_block = 0
current_var = 0

iread = 1

; output

vpixels = 600
hpixels = 800

; create a structure to hold the output information
output = {type:outputType,         $ ; flag for output medium
          hsize:float(hpixels[0]), $ ; # of pixels for graphic in x
          vsize:float(vpixels[0])} ; # of pixels for graphic in y


; plotting options
opts = [0, $ ; take the log of data
        0, $ ; take abs value of var
        0, $ ; plot max of variable
        1, $ ; draw block boundaries
        0, $ ; plot time and credits
        1, $ ; plot colorbar
        1]   ; display tick marks

; make the options the same as xflash3 defaults
blockDistLevel = 0
showProcNums = 0
options = {blocks:opts[3],   $  ; draw block boundaries
           procdist:blockDistLevel, $ ; color blocks based on processor  
           showprocnums:showProcNums, $ ; show the processor numbers on the blocks
           log:opts[0],      $  ; take the log of data
           annotate:opts[4], $  ; plot time and credits
           colorbar:opts[5], $  ; plot colorbar
           max:opts[2],      $  ; plot max of variable
           abs:opts[1],      $  ; take abs value of var
           tick:opts[6]}        ; display tick marks



; data range
dataMin = 0.1
dataMax = 0.1
auto = 1            

           
variable1 = {name:varname1,  $
            min:dataMin, $
            max:dataMax, $
            auto:1}

variable2 = {name:varname2,  $
            min:dataMin, $
            max:dataMax, $
            auto:1}


; create a structure to hold the zoom information
zoom = {plane:0, $
        xmin:-1,   $
        xmax:-1,   $
        ymin:-1,   $
        ymax:-1,   $
        zmin:-1,   $
        zmax:-1}

; a stand-in for the full widget-structure in xflash3.pro
particle_widget = {enabled:0, $
                   plot_vel:1, $
                   sym_size: 1.0, $
                   show_tag: 1, $
                   typical_velocity:10.}
            
 
problemInfo = {orientation:0, $
               name:'current default name'}

;help, data_ptr1
if (ndim EQ 2) then begin
    xplot2d_amr_diff, FILENAME1 = filename1, $
      FILENAME2 = filename2, $
      READ = iread, $
      VARIABLE_INFO1 = variable1, $
      VARIABLE_INFO2 = variable2, $
      OPTIONS = options, $
      OUTPUT = output, $
      VECTOR = vector, $
      ZOOM = zoom, $
      COLORMAP = 'Fire', $
      WIDGET_INFORMATION = info, $
      CONTOUR_OPT = contours, $
      PARTICLE_WIDGET = particle_widget, $
      PROBLEM_INFO = problemInfo, $
;      KNOWN_VARIABLES1 = varnames1, $
;      TREE_PTR = tree_ptr1, $
;      PARAMS_PTR = params_ptr1, $
;      DATA_PTR = data_ptr1, $
      PLOT_COUNT = '1', $
      DEBUG = debug, $
      NOINTERP = nointerp
;      PLOT_DATA = diffed_data


endif

	   
end
