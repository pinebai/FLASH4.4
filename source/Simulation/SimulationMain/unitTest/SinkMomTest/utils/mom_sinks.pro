pro mom_sinks

  sinkstab = dblarr(16,1000000)
  
  line = dblarr(16)
  openr, 1, 'sinks_evol.dat_cleaned'
  readf, 1, format='(A)'
  linenum = long(0)
  while not eof(1) do begin
     readf, 1, line
     sinkstab[*,linenum] = line
     linenum += 1
  endwhile
  close, 1

  sinkstab = sinkstab[*,0:linenum-1]

  id = sinkstab[0,*]
  t  = sinkstab[1,*]
  m  = sinkstab[14,*]
  x  = sinkstab[2,*]
  y  = sinkstab[3,*]
  z  = sinkstab[4,*]
  vx = sinkstab[5,*]
  vy = sinkstab[6,*]
  vz = sinkstab[7,*]
  px = vx*m
  py = vy*m
  pz = vz*m
  ; L computed for center at (0,0,0)
  lx = (y*pz-z*py) + sinkstab[11,*]
  ly = (z*px-x*pz) + sinkstab[12,*]
  lz = (x*py-y*px) + sinkstab[13,*]


  uniq_ts = t[UNIQ(t, SORT(t))]
  ntimes = n_elements(uniq_ts)

  n_tot  = lonarr(ntimes)
  m_tot  = dblarr(ntimes)
  px_tot = dblarr(ntimes)
  py_tot = dblarr(ntimes)
  pz_tot = dblarr(ntimes)
  lx_tot = dblarr(ntimes)
  ly_tot = dblarr(ntimes)
  lz_tot = dblarr(ntimes)

  for it = 0, ntimes-1 do begin

     index = where(t eq uniq_ts[it])

     id_index = id[index]
     uniq_ids = id_index[UNIQ(id_index, SORT(id_index))]
;     print, 'it= ', it, '      uniq_ids= ', uniq_ids
     if (n_elements(uniq_ids) ne n_elements(index)) then begin & print, 'uups.', uniq_ts[it], uniq_ids & stop & endif
     
     n_tot [it] = n_elements(uniq_ids)
     m_tot [it] = total(m [index],/double)
     px_tot[it] = total(px[index],/double)
     py_tot[it] = total(py[index],/double)
     pz_tot[it] = total(pz[index],/double)
     lx_tot[it] = total(lx[index],/double)
     ly_tot[it] = total(ly[index],/double)
     lz_tot[it] = total(lz[index],/double)

  endfor


  ; write file
  outfile = 'mom_sinks.dat'
  openw, 1, outfile
  printf, 1, format='(9A16)', 'time', 'n_tot', 'm_tot', 'px_tot', 'py_tot', 'pz_tot', 'lx_tot', 'ly_tot', 'lz_tot'
  for it = 0, ntimes-1 do begin
     printf, 1, format='(E16.6,I16,7E16.6)', uniq_ts[it], n_tot[it], m_tot[it], px_tot[it], py_tot[it], pz_tot[it], lx_tot[it], ly_tot[it], lz_tot[it]
  endfor
  close, 1
  print, outfile, ' written.'

  return

end
