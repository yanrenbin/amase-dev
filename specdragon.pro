pro specdragon, lambda, npixel=npixel, pixelsize=pixelsize, fcoll=fcoll, fcam=fcam, fiber=fiber, refindex=refindex,camspeed=camspeed,delta_n=delta_n,resolution=resolution,d1=d1
; lambda in Angstrom
; fcoll in mm
; fcam in mm
; pixelsize in micron
; fiber size in micron (default fiber = 50)

   lambda = float(lambda)
   if NOT keyword_set(fiber) then fiber = 50.  ; we will only use 50micron fibers.
   if NOT keyword_set(refindex) then refindex=1.5
   if NOT keyword_set(delta_n) then delta_n = 0.1
   if NOT keyword_set(d1) then d1=20
   if NOT keyword_set(resolution) then message,'Need to specify resolution keyword.'
   sampling = fiber*(float(fcam)/float(fcoll))/float(pixelsize)
   print,'Pixel sampling:', sampling, ' pixels per LSF.'
   delta_lambda = lambda/float(resolution)/float(sampling)*float(npixel)
   print,'Spectral coverage:', lambda-delta_lambda/2.,' to ',lambda+delta_lambda/2.,' A'

   vph, lambda, resolution=resolution, fcoll=fcoll, fiber=fiber, alpha=alpha, density=density,refindex=refindex

   w2 = fcam/camspeed

   w1arr = fcoll/(findgen(16)/10.+3)
   w2arr = beamwidth(w1arr,d1=d1,fcoll=fcoll,alpha=alpha,sampling=sampling,npixel=npixel,fiber=fiber,refindex=refindex,halftheta=halftheta)
   w1 = interpol(w1arr,w2arr,w2)
   collspeed= fcoll/w1

   print,'For R='+string(resolution,format='(i0.0)')+' @ '+string(lambda,format='(i0.0)')+'A, collspeed needed is : ',collspeed
   print,'Beam widens from ',w1, ' to ', w2, ' mm'

   thickness_arr = 1+findgen(191)/10.
;   delta_n = 0.1
   order = [1,2,3,4]

   angularwidth = 2*halftheta/refindex
   possible = gratingthickness(lambda,density=density/order,angularwidth=angularwidth,delta_n=delta_n,alpha=alpha,order=order) 
   if total(possible) gt 0 then begin
      print,'R = '+string(resolution,format='(i0.0)')+': Good.' 
      ind = where(possible, ct_good)
      for k=0,ct_good-1 do begin
         print,'  Order ',order[ind[k]], ' with density ',density/order[ind[k]],' l/mm'
         bestthick=optimalgrating(lambda,density=density/order[ind[k]],angularwidth=angularwidth,alpha=alpha,delta_n=delta_n,best_eta_s=best_eta_s,best_eta_p=best_eta_p,order=order[ind[k]])
	 print,'    Best thickness: ',bestthick, ' with efficiencies: ',best_eta_s,best_eta_p
      endfor
   endif else print,'R = '+string(resolution,format='(i0.0)')+': Bad.'

end
