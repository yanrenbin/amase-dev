; Output:  alpha,density
pro vph,lambda,resolution=resolution, refindex=refindex, fcoll=fcoll, fiber=fiber,alpha=alpha,density=density
; Given a required resolution at a wavlength (Angstrom), 
; and the focal length and refractive index of the substrate, 
; compute the required grating line density.

  if NOT keyword_set(resolution) then resolution = 8000.
  if NOT keyword_set(refindex) then refindex =1.5
  if NOT keyword_set(fcoll) then fcoll = 105. ; mm
  if NOT keyword_set(fiber) then fiber = 50. ; micron

;  dlambda_per_pixel = float(lambda)/float(resolution)/sampling ; Angstrom/pixel
;  lineardispersion = pixel/dlambda_per_pixel ; micron/Angstrom
;  angulardispersion = lineardispersion/(fl*1000.)*1.e7 ; radian/mm
  angulardispersion = float(resolution)/(float(lambda)*1.e-4)*float(fiber)/float(fcoll) ; radian/mm
  angulardispersion_inside = angulardispersion/refindex

  alpha = atan(float(lambda)*1.e-7*angulardispersion_inside/2.)
  density = refindex*cos(alpha)*angulardispersion_inside

;  theta1 = maxlambda-
;  density =600.+findgen(5401.)
;  sin_alpha = density*float(lambda)*1.e-7/(2*refindex)
;  k=where(sin_alpha lt 0.995)
;  alpha = asin(sin_alpha[k])
;  angdisparr= density[k]/refindex/cos(alpha)
;  finaldensity = interpol(density[k],angdisparr,angulardispersion)
;  finalalpha = interpol(alpha,angdisparr,angulardispersion)*!radeg
  print
  print,'Resolution: '+string(fix(resolution),format='(i0.0)')+' @ '+string(fix(lambda),format='(i0.0)')+'A'
;  print,'Collimator focal length: '+string(fix(fcoll),format='(i0.0)')+'mm'
;  print,'Fiber size: '+string(fix(fiber),format='(i0.0)')+'micron'
;  print,'Refractive Index: '+string(refindex,format='(f4.2)')
  print,'Required line density:',density,' l/mm'
  print,'Incidence angle within substrate:',alpha*!radeg,' deg'
  print,'Size of VPH grating assuming f/2.8 beam:',fcoll/cos(alpha)/2.8,' mm'
end

    


