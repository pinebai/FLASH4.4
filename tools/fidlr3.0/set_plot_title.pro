pro set_plot_title,TITLE=titleName,VARIABLE_INFO=variable,OPTIONS=options

; remove the setup information about creating a plot title from
; xplotNd_amr

;-------------------------------------------
; plot title setup
;--------------------------------------------
; Set the long version of the variable name, for the plot title
titleName = variable.name

; override the short name for the standard cases
case variable.name of

    'dens': titleName = 'Density (g/cm!E3!N)'
    'temp': titleName = 'Temperature (K)'
    'gamc': titleName = 'Gamma C'
    'game': titleName = 'Gamma E'
    'velx': titleName = 'X Velocity (cm/s)'
    'vely': titleName = 'Y Velocity (cm/s)'
    'velz': titleName = 'Z Velocity (cm/s)'
    'pres': titleName = 'Pressure (erg/cm!E3!N)'
    'ener': titleName = 'Energy (ergs/g)'
    'enuc': titleName = 'Nuclear Energy Generation Rate (ergs/g/s)'

    'tot_vel':   titleName = 'Total Velocity (cm/s)'
    'int_ener':  titleName = 'Specific Internal Energy (erg/g)'
    'ekin/eint': titleName = 'Kinetic Energy / Internal Energy'
    'snd_spd':   titleName = 'Sound Speed (cm/s)'
    'mach':      titleName = 'Mach Number'

    else:

endcase

; deal with the xplot2d_amr_diff case, where all these options might
;  not exist
IF (N_Elements(options.abs) NE 0) AND (options.abs) THEN titleName = 'Abs ' + titleName
if (N_Elements(options.log) NE 0) AND (options.log) THEN  titleName = 'Log10 ' + titleName
if (N_Elements(options.max) NE 0) AND (options.max) THEN  titleName = 'Max ' + titleName


if N_Elements(options.procdist NE 0) then BEGIN
    IF (options.procdist) THEN BEGIN 
        titleName = 'Proc Distribution'
        options.colorbar = 0
    ENDIF 
ENDIF

; end of plot title setup

END 
