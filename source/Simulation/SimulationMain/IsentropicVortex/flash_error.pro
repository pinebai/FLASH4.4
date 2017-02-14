;
; This file contains the following idl procedures and functions:
;   pro      numerror
;   pro      error_table_output
;   function convergenceRate2d
;   function convergenceRate1d
;   pro      flash_error
; The one you want to use is (probably) flash_error, which calls
; the others. numerror calls loaddata, a procedure in fidlr2. 
; Before running numerror or flash_error, you need to compile start.pro
; or start_linux.pro to access loaddata.
;
; More details are given in the comments for flash_error, below.
;


;  
;   pro numerror, exactfile, testfile,
;                 VAR=var,
;                 ERROR=error,
;                 texact, ttest,
;                 gridsizes,
;                 enorms
;
; Find the error in a numerical solution. Actually, Find the
; difference between two numerical solutions, one of which
; is assumed to be the exact solution.
;
; Assumptions:
;  (1) equispaced cartesian grid
;  (2) your exact and test solutions really should match up -
;        if the test solution time is slightly off or the grids
;        don't overlap exactly, all bets are off.
;
;       Input
;         exactfile -- name of the FLASH file to read in, which contains
;            the "exact solution"
;         testfile  -- name of the FLASH file to read in, which contains
;            the solution to be tested
;         var       -- the variable to be used to compute the error
;       Output
;         error     -- array of error at each grid point
;         texact    -- time of exact solution
;         ttest     -- time of test solution
;         gridsizes -- array holding number of dimensions, nx, and ny
;                        for the exact and test solutions
;             gridsizes[0]  ndim
;             gridsizes[1]  nx
;             gridsizes[2]  ny
;             gridsizes[3]  ?
;             gridsizes[4]  nx*ny
;         enorms    -- array holding computed error norms
;             enorms[0]  infinity-norm
;             enorms[1]  1-norm
;             enorms[2]  2-norm


pro numerror, exactfile, testfile, VAR=var, ERROR=error, $
              texact, ttest, gridsizes, enorms

    if (not(keyword_set(var)))  then var = 'dens'

; ----- Read the "exact" solution -----
;    loadfile, exactfile, exactvar, VARNAM=var,   $
;              XCOORDS=xexact, YCOORDS=yexact, ZCOORDS=zexact, FILETIME=texact
    exactvar = loaddata(exactfile, var,                                   $
               XCOORDS=xexact, YCOORDS=yexact, ZCOORDS=zexact, TIME=texact)
    gridsizes = size(exactvar)
    
; ----- Read the test solution -----
;    loadfile, testfile,   testvar, VARNAM=var,   $
;              XCOORDS=xtest,  YCOORDS=ytest, ZCOORDS=ztest,  FILETIME=ttest
    testvar = loaddata(testfile, var,                                   $
              XCOORDS=xtest,  YCOORDS=ytest, ZCOORDS=ztest,  TIME=ttest)
    testsizes = size(testvar)

; ----- get information about solutions  -----
; ----- including some error checking... -----
    ndim = gridsizes[0]
    case ndim of
      1: if (testsizes[0] NE 1) then begin 
            print, "numerror error: the dimensions of ", testfile
            print, "  must match the dimensions of ", exactfile
         endif
      2: if (testsizes[0] NE 2) then begin 
            print, "numerror error: the dimensions of ", testfile
            print, "  must match the dimensions of ", exactfile
         endif
      3: if (testsizes[0] NE 3) then begin 
            print, "numerror error: the dimensions of ", testfile
            print, "  must match the dimensions of ", exactfile
         endif
      else: print, "numerror error: can't handle ", ndim, " dimensions!"
    endcase
;
    dsizes = gridsizes - testsizes
    for j=1, ndim do begin
      if (dsizes[j] NE 0) then begin
         print, "dimension ",j,": ",exactfile, " and ", testfile
         print, "   are not the same size!"
      endif
    endfor

;
; ----- now compute the error norms -----
;   This is the difference between the two solutions.
    error = testvar - exactvar

;   Assume an equispaced grid here.
    case ndim of
      1: dvol =  xexact[2]-xexact[1]
      2: dvol = (xexact[2]-xexact[1])*  $
                (yexact[2]-yexact[1])
      3: dvol = (xexact[2]-xexact[1])*  $
                (yexact[2]-yexact[1])*  $
                (zexact[2]-zexact[1])
      else: print, "numerror error: can't handle ", ndim, " dimensions!"
    endcase

;   Calculation of error norms.
    einf = max( abs(error) )*dvol
    e1 = total( abs(error), /double)*dvol
    e2 = sqrt( total(error^2, /double)*dvol )
    enorms = dindgen(3)
    enorms[0] = einf
    enorms[1] = e1
    enorms[2] = e2

end

;
;   In 2d, estimate h_coarse/h_fine by taking the square root of the
;     coarse array size/fine array size and inverting.

function convergenceRate2d, coarse_sizes,  fine_sizes,  $
                            coarse_errors, fine_errors


    size_ratio = 1.0d0/sqrt( double(coarse_sizes[4])/double(fine_sizes[4]) )

    rate = alog10(coarse_errors/fine_errors)/alog10(size_ratio)

    return, rate

end

;
;   In 1d, estimate h_coarse/h_fine by taking the
;     coarse array size/fine array size and inverting.

function convergenceRate1d, coarse_sizes,  fine_sizes,  $
                            coarse_errors, fine_errors


    size_ratio = 1.0d0/( double(coarse_sizes[3])/double(fine_sizes[3]) )

    rate = alog10(coarse_errors/fine_errors)/alog10(size_ratio)

    return, rate

end

;
; pro error_table_output, output_file, ngrids,
;                         dimensions, var,
;                         efiles, tfiles,
;                         gsizes, enorms, rates
;

pro error_table_output, output_file, ngrids,     $
                        dimensions, var,         $
                        efiles, tfiles,          $
                        gsizes, enorms, rates

;
; ----- format statements -----
;
; assume that filenames are no more than 34 characters long. two
;   filenames will fit on an 72 character line.

  filelist_format = '(2(2x,A34))'
  case dimensions of
    1: begin 
         line1_format = '(i4,3(1x,e12.5,2x, 7x ))'
         lines_format = '(i4,3(1x,e12.5,2x,f7.5))'
       end
    2: begin 
         line1_format = '(i4,1x,i4,3(1x,e12.5,2x, 7x ))'
         lines_format = '(i4,1x,i4,3(1x,e12.5,2x,f7.5))'
       end
    else: print, 'error_table_output error: unrecognized dimensions!"
  endcase

;
; -----        open the output file       -----
; ----- write information about the tests -----
;
  openw, output_lun, output_file, /get_lun

  printf, output_lun, '----------------------------------------'

  case dimensions of
    1: printf, output_lun, '1d error comparison on ', var
    2: printf, output_lun, '2d error comparison on ', var
    else: print, 'error_table_output error: unrecognized dimensions!"
  endcase

  printf, output_lun, ' '

  printf, output_lun, 'Exact solutions:                      ',  $
                      'Test solutions:                       ',  $
                       format=filelist_format

  for j=0, ngrids-1 do begin
    printf, output_lun, efiles[j],                     $
                        tfiles[j],                     $ 
                        format=filelist_format
  endfor

  printf, output_lun, ' '

;
; ----- now write the table -----
;

  case dimensions of
    1: begin
;        Table header line
         printf, output_lun, '  nx ',                        $
                             '    E_inf   ', '   Rate   ',   $
                             '     E_1    ', '   Rate   ',   $
                             '     E_2    ', '   Rate '
;        First line
         printf, output_lun, gsizes[1,0],                      $
                             enorms[0,0],                      $
                             enorms[1,0],                      $
                             enorms[2,0],                      $
                             format=line1_format
;        Rest of lines
         for j=1, ngrids-1 do begin
           printf, output_lun, gsizes[1,j],                $
                               enorms[0,j], rates[0,j-1],  $
                               enorms[1,j], rates[1,j-1],  $
                               enorms[2,j], rates[2,j-1],  $
                               format=lines_format
         endfor
       end

    2: begin
;        Table header line
         printf, output_lun, '  nx ', '  ny ',               $
                             '    E_inf   ', '   Rate   ',   $
                             '     E_1    ', '   Rate   ',   $
                             '     E_2    ', '   Rate '
;        First line
         printf, output_lun, gsizes[1,0],   gsizes[2,0],   $
                             enorms[0,0],                  $
                             enorms[1,0],                  $
                             enorms[2,0],                  $
                             format=line1_format
;        Rest of lines
         for j=1, ngrids-1 do begin
           printf, output_lun, gsizes[1,j], gsizes[2,j],   $
                               enorms[0,j], rates[0,j-1],  $
                               enorms[1,j], rates[1,j-1],  $
                               enorms[2,j], rates[2,j-1],  $
                               format=lines_format
         endfor
       end

    else: print, 'error_table_output error: unrecognized dimensions!"

  endcase

  printf, output_lun, '----------------------------------------'

;
; ----- all done -----
;
  free_lun, output_lun

end

;
;  pro flash_error, INPUT_FILE=input_file, OUTPUT_FILE=output_file,
;    GRID_SIZES=run_gsizes, ENORMS=run_enorms, RATES=run_rates
;
; A procedure to compute convergence rates for a series of FLASH runs.
;
; The input and output are mainly handled through files, so two parameters
;   to the procedure are the file names. In addition, some results can be
;   passed out of the procedure for more analysis in IDL (e.g. plotting)
;
;
; These are inputs to the procedure:
; input_file  == the name of the input file; if not given, defaults
;                  to flash_error.in
;                The input file format is as follows:
;                  int ngrids   == the number of grids to compute the error on
;                  string efile == the file(s) holding the exact solution,
;                                    one file per line (there should be ngrids
;                                    files)
;                  string tfile == the file(s) holding the test solution,
;                                    one file per line (there should be ngrids
;                                    files)
;                  string var   == the variable to be tested. Use the FLASH
;                                    four letter name, eg dens, velx, pres, etc.
;
; output_file == the name of the output file; if not given, defaults
;                  to flash_error.out
;
; These are procedure outputs to the calling program/idl session:
; run_gsizes  == an array containing data about the solutions tested
;
; run_enorms  == an array containing the computed error norms
;
; run_rates   == an array continaing the computed convergence rates
;
;
; This procedure uses the procedures numerror and error_table_output and
; the function convergenceRate2d, all included below.
;

pro flash_error, INPUT_FILE=input_file, OUTPUT_FILE=output_file, $
  GRID_SIZES=run_gsizes, ENORMS=run_enorms, RATES=run_rates

  if (not(keyword_set(input_file)))   then input_file  = 'flash_error.in'
  if (not(keyword_set(output_file)))  then output_file = 'flash_error.out'

  ngrids = 0
  testvar = ''

;
; Read the input file.
; Read the number of grids.
; Make some declarations based on ngrids. 
; Read the the file names of the exact solutions on those grids,
;          the file names of the test solutions,
;          the variable to be tested.

  openr, input_lun, input_file, /get_lun

  readf, input_lun, ngrids

  print, 'Number of grids: ',ngrids

  exactfile = ''
  testfile = ''
  run_efile = strarr(ngrids)
  run_tfile = strarr(ngrids)

  for j=0, ngrids-1 do begin
    readf, input_lun, exactfile
    run_efile[j] = exactfile
    print, 'Grid ',j, ', Exact solution file ', run_efile[j]
  endfor
  for j=0, ngrids-1 do begin
    readf, input_lun, testfile
    run_tfile[j] = testfile
    print, 'Grid ',j, ', Test solution file ', run_tfile[j]
  endfor

  readf, input_lun, testvar

  print, 'Test variable: ', testvar

  free_lun, input_lun

;
; More declarations:
;

  run_ttest  = dblarr(ngrids)
  run_gsizes = ulon64Arr(5,ngrids) ; 5 : size of the output of idl fcn "sizes".
  run_enorms = dblarr(3,ngrids)    ; 3 : the number of norms (inf, 1, and
  run_rates  = dblarr(3,ngrids-1)  ;    2-norm). then, one rate for each norm.

;
; Now the work:
;

  for j=0,ngrids-1 do begin  ; j is each grid size 

    exactfile = run_efile[j]
     testfile = run_tfile[j]

    numerror, exactfile, testfile,VAR=testvar,         $
              t_exact, t_test, gridsizes, enorms

    run_ttest[j] = t_test
    run_gsizes[*,j] = gridsizes
    run_enorms[*,j] = enorms

  endfor

  ndim = run_gsizes[0,0]
  for k=0,ngrids-2 do begin  ;  k is "between" each grid size

    coarse_sizes  = run_gsizes[*,k  ]
    fine_sizes    = run_gsizes[*,k+1]
    coarse_enorms = run_enorms[*,k  ]
    fine_enorms   = run_enorms[*,k+1]

    case ndim of
      1: run_rates[*,k] = convergenceRate1d( coarse_sizes,  fine_sizes,   $
                                             coarse_enorms, fine_enorms )
      2: run_rates[*,k] = convergenceRate2d( coarse_sizes,  fine_sizes,   $
                                             coarse_enorms, fine_enorms )
      default: print, "Sorry, didn't implement convergenceRate for ",      $
           ndim, " dimensions yet."
    endcase

  endfor

;
; Write output to a file.
;

  ndim = run_gsizes[0,0]
  error_table_output, output_file, ngrids,                  $
                      ndim, testvar, run_efile, run_tfile,  $
                      run_gsizes, run_enorms, run_rates

end

