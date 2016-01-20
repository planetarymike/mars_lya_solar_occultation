pro sav2cppinput, fname, outfilename=outfilename
  if n_elements(outfilename) eq 0 then outfilename='test_out.dat'
  
  ;read in the data from the save file
  restore, fname

  ;;based on an email from Ed, "EUVM H corona occultation study", 14 Aug 2015

  ;; The variables are:
  ;; h_tangent is the height of the LOS vector above the limb in km.
  ;; h_sc is the MAVEN height.
  ;; SZA is the solar zenith angle in degrees.
  ;; i_cnts is the lyman alpha irradiance in w/m^2
  ;; And i_max is the maximum signal.

  nscalt = n_elements(h_sc)
  if n_elements(h_tangent) NE nscalt or $
     n_elements(SZA)       NE nscalt or $
     n_elements(i_cnts)    NE nscalt    $
  then begin
     print, "number of elements in sav file is not shared across components! Aborting"
     print,  thisisonlyheretocauseanerror[404]
  endif

  ;stop
  
  openw, 1, outfilename
  printf, 1, nscalt, FORMAT = '(I-)'
  for i=0,nscalt-1 do begin
     printf, 1, $
             [h_sc[i], h_tangent[i], -SZA[i]], $
             FORMAT = '(2(F14.2," "),(F11.5," "))'
  endfor
  close,1
  
end