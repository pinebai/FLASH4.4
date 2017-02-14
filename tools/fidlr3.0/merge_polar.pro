; merge a FLASH AMR dataset in polar coordinates onto a uniform grid

; 7-10-03
;
; Michael Zingale
; University of California, Santa Cruz
;
; *** not for redistribution ***
;
;
; the basic steps are the following:
;
; 1. find all the blocks were a going to consider, but picking the
;    highest resolution blocks at a given coordinate where the area of
;    a zone is >~ the area of a pixel in our uniform data array.
;
; 2. consider each uniform data zone, and subdivide it.  
;
; 3. convert the center of each subzone to r, theta, and find the data
;    value it corresponds to -- if there is more than one data/polar
;    zone it corresponds to, then use the finer one.
;
; 4. average all the subzones and store -- be careful here -- we may
;    be at the boundary of the domain, and not all subzones may have
;    data.  This needs to be taken into account when doing the
;    averaging.  We will only fill a uniform data cell if >= 1/2 of
;    the subzones have data.
;

function merge_polar, temp_arr, $
                      TREE = tree, PARAMETERS = params, $
                      XMERGE = xout, YMERGE = yout, $
                      SAMPLE = sample, $
                      RRANGE = rrange, $
                      TRANGE = trange

NPTS_MAX = 512
NPTS_MIN = 256

; find the minimum and maximum refinement levels of the good data
index_good = where(tree[*].nodetype EQ 1)

lmax = max(tree[index_good].lrefine)
lmin = min(tree[index_good].lrefine)

lwant = lmax

theta_min_global = min(tree[*].bndBox[0,1])
theta_max_global = max(tree[*].bndBox[1,1])

if n_elements(trange) EQ 0 then begin
    theta_min = theta_min_global
    theta_max = theta_max_global
endif else begin
    theta_min = trange[0] > theta_min_global
    theta_max = trange[1] < theta_max_global
endelse

dtheta = theta_max - theta_min

; for now, assume that we start at r = 0
r_min_global = min(tree[*].bndBox[0,0])
r_max_global = max(tree[*].bndBox[1,0])

if n_elements(rrange) EQ 0 then begin
    r_min = r_min_global
    r_max = r_max_global
endif else begin
    r_min = rrange[0] > r_min_global
    r_max = rrange[1] < r_max_global
endelse

 
; figure out how big of a uniformly gridded data array we need.  For 
; simplicity, we will always put the theta = 0 vertically.  Depending
; on how many quadrants were span, our x and y extrema change

if (dtheta LT !pi/2. - 1.e-4) then begin

; we are less than one quadrant

    if (theta_max LE !pi/2) then begin
        xmin = 0.d0
        ymin = 0.d0

        xmax = r_max*sin(theta_max)
        ymax = r_max*cos(theta_min)
        
    endif else if (theta_max LT !pi) then begin
        xmin = 0.d0
        ymax = 0.d0

        xmax = r_max*sin(theta_min)
        ymin = r_max*cos(theta_max)

    endif else begin
        print, 'ERROR: angle range not yet suported'
        return, -1        
    endelse

endif else begin

; we span more than one quadrant

    if ( theta_max GT !pi ) then begin

       xmin = -r_max
       xmax =  r_max

       ymin = -r_max
       ymax =  r_max

    endif else begin

       xmin = 0.d0
       xmax = r_max

       ymin = r_max*cos(theta_max)
       ymax = r_max*cos(theta_min)

    endelse

endelse

; we now have the coordinates of the Cartesian array that we are going
; to uniformly grid our data into.  So let's go ahead and allocate it.
; Eventually we will have an option to do double or float.

; dx and dy are going to be set by the radial spacing.  Eventually, we
; should perhaps consider the theta spacing too.
dr = 0.33333d0*(r_max_global - r_min_global)/(params.ntopx*params.nxb*2^(lwant-1))

if ( (long((r_max - r_min)/dr) > 1) LT NPTS_MIN) then begin

; don't let there be less than NPTS_MIN radial zones on the plot
    dr = (r_max - r_min)/double(NPTS_MIN)

endif else if ( (long((r_max - r_min)/dr) > 1) GT NPTS_MAX) then begin
    
; also don't let there be more than NPTS_MAX zones on the plot -- this is 
; also somewhat arbitrary.  For a graphics display, this is probably
; bigger than the plot will be (since we are talking just the plot
; size, not counting the axis, labels, ...  For a postscript device,
; we probably need something bigger.
    dr = (r_max - r_min)/double(NPTS_MAX)

endif


nx = long((xmax - xmin)/dr) > 1
ny = long((ymax - ymin)/dr) > 1

cartData = dblarr(nx,ny)

; fill the array with some knowingly wrong data
cartData[*,*] = -1.e33

; and also compute the vectors of x and y coordinates for the
; Cartesian gridded data.
dx = dr
xl = dindgen(nx)*dx + xmin
xr = (dindgen(nx) + 1.0)*dx + xmin
x = 0.5*(xl + xr)

dy = dr
yl = dindgen(ny)*dy + ymin
yr = (dindgen(ny) + 1.0)*dy + ymin
y = 0.5*(yl + yr)

;debug
;help, xl, xr
;help, yl, yr
;print, xl[0],xr[nx-1]
;print, yl[0],yr[ny-1]

; ok... now the easy stuff is done.  Now comes the hard stuff.  

; figure out which blocks it is we are going to plot.  We want to make
; sure the area of the block is not too much smaller than the area of
; the Cartesian zone.  We will actually do this on a block by block 
; basis.
cartArea = dr*dr

rl = dblarr(params.nxb, params.totBlocks)
rr = dblarr(params.nxb, params.totBlocks)

tl = dblarr(params.nyb, params.totBlocks)
tr = dblarr(params.nyb, params.totBlocks)

;debug
;help, rl, rr
;help, tl, tr

; we need to know the extrema of all the polar zones.  This probably
; won't scale all that well, but we'll deal with that later.
for i = 0l, params.totBlocks-1 do begin


    rlcoords = (tree[i].bndBox[1,0] - tree[i].bndBox[0,0])* $
      ((dindgen(params.nxb))/double(params.nxb)) + tree[i].bndBox[0,0]

    rrcoords = (tree[i].bndBox[1,0] - tree[i].bndBox[0,0])* $
      ((dindgen(params.nxb) + 1.0)/double(params.nxb)) + tree[i].bndBox[0,0]
    

    tlcoords = (tree[i].bndBox[1,1] - tree[i].bndBox[0,1])* $
      ((dindgen(params.nyb))/double(params.nyb)) + tree[i].bndBox[0,1]

    trcoords = (tree[i].bndBox[1,1] - tree[i].bndBox[0,1])* $
      ((dindgen(params.nyb) + 1.0)/double(params.nyb)) + tree[i].bndBox[0,1]



    rl[*,i] = rlcoords
    rr[*,i] = rrcoords

    tl[*,i] = tlcoords
    tr[*,i] = trcoords

;debug
;   print, rl[0,i],rr[7,i],tl[0,i],tr[7,i]

endfor

; needs to be modified

blkArea = !pi*(tree[*].bndBox[1,0] - tree[*].bndBox[0,0])* $
              (tree[*].bndBox[1,0] + tree[*].bndBox[0,0])* $
              (tree[*].bndBox[1,1] - tree[*].bndBox[0,1])

; now figure out which blocks we are going to actually try to use to
; fill the cartesian gridded data.  The trick is that if the block is
; too small, then it is probably not worth plotting.  When the data
; was stored, it was restricted up the tree in a conservative fashion,
; and this is is more accurate then the sampling that we are going to
; be doing here.

; we need to make sure that we don't exclude any blocks, so the list
; of good blocks must include all base level blocks

; if you try to play smart and pick blocks based on their area, the
; icart where line will occupy about 75% of your runtime.


; ** OPTIMIZATION ** 
; getting this routine to run fast is all about getting the list of 
; good blocks as small as possible while not throwing away the 
; important ones.

;igood = where(blkArea GE 0.25*params.nxb*params.nyb*cartArea OR $
;                     tree[*].lrefine EQ 1)

;igood = where(tree[*].nodeType EQ 1)

; remove any blocks that are outside the zoom region

igood1 = where(tree[*].bndBox[1,0] GE r_min AND $
              tree[*].bndBox[0,0] LE r_max*sqrt(2.d0) AND $
              tree[*].bndBox[1,1] GE theta_min AND $
              tree[*].bndBox[0,1] LE theta_max AND $
              tree[*].nodeType EQ 1)

igood2 = where(tree[*].bndBox[1,0] GE r_min AND $
              tree[*].bndBox[0,0] LE r_max*sqrt(2.d0) AND $
              tree[*].bndBox[1,1] GE theta_min AND $
              tree[*].bndBox[0,1] LE theta_max AND $
              (tree[*].lrefine EQ 1 OR $
               blkArea GE params.nxb*params.nyb*cartArea))


if ((size(igood1))[1] GT (size(igood2))[1]) then begin
    igood = igood2
endif else begin
    igood = igood1
endelse

;debug
;help, igood
;print, igood

; now loop over all of the *Cartesian* cells and start a filling them.
; We are going to be doing subsampling to build up a (semi-)accurate
; average.  Eventually, this subzoning size will be variable.


; make a two dimensional x, y space
xl_ij = xl # replicate(1.d0, ny)
xr_ij = xr # replicate(1.d0, ny)

yl_ij = replicate(1.d0, nx) # yl
yr_ij = replicate(1.d0, nx) # yr

r_cart = dblarr(4,nx,ny)

r_cart[0,*,*] = xl_ij^2 + yl_ij^2
r_cart[1,*,*] = xl_ij^2 + yr_ij^2
r_cart[2,*,*] = xr_ij^2 + yl_ij^2
r_cart[3,*,*] = xr_ij^2 + yr_ij^2

t_cart = dblarr(4,nx,ny)

t_cart[0,*,*] = atan(xl_ij, yl_ij)
t_cart[1,*,*] = atan(xl_ij, yr_ij)
t_cart[2,*,*] = atan(xr_ij, yl_ij)
t_cart[3,*,*] = atan(xr_ij, yr_ij)

;debug
;print, nx, ny
;help, t_cart

for j = 0l, ny-1 do begin
   for i = 0l, nx-1 do begin
      if  (t_cart[0,i,j] LT 0) then begin
           t_cart[0,i,j] = 2.d0 * !pi + t_cart[0,i,j]
      endif
      if  (t_cart[1,i,j] LT 0) then begin
           t_cart[1,i,j] = 2.d0 * !pi + t_cart[1,i,j]
      endif
      if  (t_cart[2,i,j] LT 0) then begin
           t_cart[2,i,j] = 2.d0 * !pi + t_cart[2,i,j]
      endif
      if  (t_cart[3,i,j] LT 0) then begin
           t_cart[3,i,j] = 2.d0 * !pi + t_cart[3,i,j]
      endif
   endfor
endfor

for j = 0l, ny-1 do begin


    r_cart_min = sqrt(min(r_cart[*,*,j]))
    r_cart_max = sqrt(max(r_cart[*,*,j]))
    
    t_cart_min = min(t_cart[*,*,j])
    t_cart_max = max(t_cart[*,*,j])
        
    icart_column = where(tree[igood].bndBox[1,0] GT r_cart_min AND $
                         tree[igood].bndBox[0,0] LT r_cart_max AND $
                         tree[igood].bndBox[1,1] GT t_cart_min AND $
                         tree[igood].bndBox[0,1] LT t_cart_max)

    cart_column_blocks = igood[icart_column]

    for i = 0l, nx-1 do begin
        

; start by figuring out if the current Cartesian zone is actually in
; the r, theta domain -- it may not be, what with the circles and
; everything.

        r_cart_min = sqrt(min(r_cart[*,i,j]))
        r_cart_max = sqrt(max(r_cart[*,i,j]))

        t_cart_min = min(t_cart[*,i,j])
        t_cart_max = max(t_cart[*,i,j])


        icart = where(tree[cart_column_blocks].bndBox[1,0] GT r_cart_min AND $
                      tree[cart_column_blocks].bndBox[0,0] LT r_cart_max AND $
                      tree[cart_column_blocks].bndBox[1,1] GT t_cart_min AND $
                      tree[cart_column_blocks].bndBox[0,1] LT t_cart_max)


        if (icart[0] NE -1) then begin

; we need the block numbers wrt to the original numbers        
            cart_blocks = cart_column_blocks[icart]


; the subzone size depends on how coarse the blocks are
            lrefine_sub_max = max(tree[cart_blocks].lrefine)

            nsub = 2 - (lmax - lrefine_sub_max) > 1

            dx_sub = dx/nsub
            dy_sub = dy/nsub
        
            data_sub = dblarr(nsub,nsub)

; we should actually do this check before the subsampling eventually
            data_sub[*,*] = 0.d0
            numSub = 0
            vol = 0.d0
        
            xc_sub = (dindgen(nsub)+0.5)*dx_sub + xl[i]
            yc_sub = (dindgen(nsub)+0.5)*dy_sub + yl[j]
        
            rc_sub = sqrt(xc_sub^2 + yc_sub^2)
            tc_sub = atan(xc_sub,yc_sub)

            for jj = 0, nsub-1 do begin
                if (xc_sub[jj] LT 0) then begin
                   tc_sub[jj] = 2.d0 * !pi + tc_sub[jj]
                endif
            endfor
        
; we need to do a weighted average by the volume of a Cartesian cell.
; Since of r, theta spherical coordinates is axisymmetric, we compute
; the volume of a subzone rotated about the vertical axis.  The
; vol_sub will only depend on x position, so it is a 1-d array
            vol_sub = 2.d0*!pi*dy_sub*xc_sub*dx_sub


            for jj = 0, nsub-1 do begin
                for ii = 0, nsub-1 do begin

; this part can probably be done in a vector fashion, up until the
; match_blocks part

                    imatch = where(rc_sub[ii] GE double(tree[cart_blocks].bndBox[0,0]) AND $
                                   rc_sub[ii] LT double(tree[cart_blocks].bndBox[1,0]) AND $
                                   tc_sub[jj] GE double(tree[cart_blocks].bndBox[0,1]) AND $
                                   tc_sub[jj] LT double(tree[cart_blocks].bndBox[1,1]))


                    if (imatch[0] NE -1) then begin

                        match_blocks = cart_blocks[imatch]
;                        match_blocks = imatch


; take only the block that is at the highest level of refinement out
; of these
                        refine_gb = tree[match_blocks].lrefine
                        
; sometimes, for some reason, there are two blocks at the same
; refinement level that ``contain'' the subzone, due to roundoff
; error.  Take the first one here.  Yes, this is a hack.
                        dataBlock = (match_blocks[where(refine_gb EQ max(refine_gb))])[0]


; ok, know we know the block which contains the current subzone, so
; figure out the cell
                        zone_i = where(rc_sub[ii] GE rl[*,dataBlock] AND $
                                       rc_sub[ii] LT rr[*,dataBlock])
                        
                        zone_j = where(tc_sub[jj] GE tl[*,dataBlock] AND $
                                       tc_sub[jj] LT tr[*,dataBlock])
                        
                        data_sub[ii,jj] = temp_arr[dataBlock,zone_i,zone_j]*vol_sub[ii]
                        vol = vol + vol_sub[ii]

                    endif
                    
                endfor
            endfor

; now do the averaging
            if (vol NE 0.0) then begin
                cartData[i,j] = total(data_sub)/vol                
            endif else begin
                cartData[i,j] = -1.e33
            endelse

        endif

    endfor
endfor

;debug
;help, x, y

xout = x
yout = y

return, cartData

end
