pro multiraster_amr, dirname 

; this routine produces a multiple raster plots of an amr object. 

; read in header
header = read_amrheader (dirname)

; get fields
logrho  = get_amrcomponent (dirname, 'log_den')

v_x  = get_amrcomponent (dirname, 'xvel')

v_y  = get_amrcomponent (dirname, 'yvel')

v_z  = get_amrcomponent (dirname, 'zvel')

rad    = get_amrcomponent (dirname, 'rad')

ieng = get_amrcomponent (dirname, 'ieng')

phi  = get_amrcomponent (dirname, header.quantities [8])

; Construct slices through equator. Output files to postscript format.

loadct, 0

; SET_PLOT, 'PS'

; DEVICE, FILENAME = 'multiraster.ps'

!P.MULTI = [0, 2, 3]

!X.STYLE = 1
!Y.STYLE = 1

raster_amr, logrho, 0, 0, /isotropic, xtitle = 'X Axis', ytitle = 'Y Axis', $
                      title = 'Equatorial Slice Log Density'

raster_amr, v_x, 0, 0, /isotropic, xtitle = 'X Axis', ytitle = 'Y Axis', $
                      title = 'Equatorial Slice X-Velocity' 

raster_amr, v_y, 0, 0, /isotropic, xtitle = 'X Axis', ytitle = 'Y Axis', $
                      title = 'Equatorial Slice Y-Velocity'

raster_amr, v_z, 0, 0, /isotropic, xtitle = 'X Axis', ytitle = 'Y Axis', $
                      title = 'Equatorial Slice Z-Velocity' 

raster_amr, rad, 0, 0, /isotropic, xtitle = 'X Axis', ytitle = 'Y Axis', $
                      title = 'Equatorial Slice Internal Energy'

raster_amr, phi, 0, 0, /isotropic, xtitle = 'X Axis', ytitle = 'Y Axis', $
                      title = 'Equatorial Slice Potential'
		      
		     
return
end



