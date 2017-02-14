PRO set_number_plots,PLOT_COUNT=plot_count, $
                     NX=nx,NY=ny,  $
                     MULTI_PLOTS=multi_plots_per_page, $
                     DEBUG = debug
;  
; sets the NX/NY variables for number of plots per page
;

case plot_count of 
    '1': begin
        NX = 1
        NY = 1
    end
    '1x2': begin
        NX=1
        NY=2
    end
    '2x1': begin
        NX=2
        NY=1
    end
    '2x2': begin
        NX = 2
        NY = 2
    end
    '2x3': begin
        NX = 2
        NY = 3
    end
    '3x2': begin
        NX = 3
        NY = 2
    end
    '3x3': begin
        NX = 3
        NY = 3
    end
ENDCASE

if (NX GT 1 OR NY GT 1) then multi_plots_per_page=1 else multi_plots_per_page=0

IF Keyword_Set(debug) THEN print, 'multi_plots_per_page = ', multi_plots_per_page


END 
