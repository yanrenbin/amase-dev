function optimalgrating,lambda,density=density,angularwidth=angularwidth,alpha=alpha,best_eta_s=best_eta_s,best_eta_p=best_eta_p,delta_n=delta_n,order=order
     spacing = 1000/density 
     if NOT keyword_set(delta_n) then delta_n = 0.1
     thickness_upper = 2*spacing/angularwidth/order*1.6
;     print,'Upper limit on thickness:',thickness_upper
     thickness_lower = asin(sqrt(0.5))/!pi/delta_n/cos(2*alpha)*(lambda*1.e-4*cos(alpha))
;     print,'Lower limit on thickness:',thickness_lower
     npix = 101
     darr = findgen(npix)/(npix-1.)*(thickness_upper-thickness_lower)+thickness_lower
     eta_s = sin(!pi*delta_n*darr/(lambda*1.e-4*cos(alpha)))^2
     eta_p = sin(!pi*delta_n*darr/(lambda*1.e-4*cos(alpha))*cos(2*alpha))^2
     plot,eta_s,eta_p,ps=1
     tmp = max(eta_s*eta_p,ind)
     deta = eta_s - eta_p
     ddeta = deta[1:npix-1]-(shift(deta,1))[1:npix-1]
     lastpoint = 0
     for i=0,npix-3 do begin
         if ddeta[i]*ddeta[i+1] le 0 or i eq npix-3 then begin
	     thispoint = i+1
	     d_option =interpol(darr[lastpoint:thispoint],deta[lastpoint:thispoint],0.0)
	     if d_option lt thickness_lower or d_option gt thickness_upper then continue
	     if lastpoint eq 0 then all_d_option= d_option else all_d_option=[all_d_option,d_option]
	     lastpoint = thispoint
	 endif
     endfor
     if n_elements(all_d_option) gt 0 then all_d_option = [darr[0],all_d_option,darr[npix-1]] else all_d_option=[darr[0],darr[npix-1]]
     eta_s_option = interpol(eta_s,darr,all_d_option)
     eta_p_option = interpol(eta_p,darr,all_d_option)
     ss = reverse(sort(eta_s_option*eta_p_option))
     best_eta_p = eta_p_option[ss[0]]
     best_eta_s = eta_s_option[ss[0]]
     t = where(eta_s_option[ss]*eta_p_option[ss] gt 0.25,ct_good)
     if ct_good gt 0 then  print,'All good options:',all_d_option[ss[t]],eta_s_option[ss[t]]
     print,'Half of Angular bandwidth produced:',spacing/all_d_option[ss[0]]/order*!radeg,' deg'
     return,all_d_option[ss[0]]
end
