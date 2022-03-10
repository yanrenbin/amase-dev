function gratingthickness,lambda,density=density,angularwidth=angularwidth,alpha=alpha,delta_n=delta_n,order=order
     spacing = 1000/density 
     if NOT keyword_set(delta_n) then delta_n = 0.1
     thickness_upper = spacing/(angularwidth/2.)/order*1.6 ; The division by is because of delta_beta = 2*delta_theta. 
;     The muliplication with 1.6 is because the half-efficiency window is 1.6 times wider than the approximate formula when nu = 1.5*pi, which is the regime we are exploring.
     print,'Trying density: ',density
     print,'Upper limit on thickness:',thickness_upper
     eta_s = sin(!pi*delta_n*thickness_upper/(lambda*1.e-4*cos(alpha)))^2
     eta_p = sin(!pi*delta_n*thickness_upper/(lambda*1.e-4*cos(alpha))*cos(2*alpha))^2
     thickness_lower = abs(asin(sqrt(0.5))/!pi/delta_n/cos(2*alpha)*(lambda*1.e-4*cos(alpha)))
     print,'Lower limit on thickness:',thickness_lower
     possible = thickness_upper gt thickness_lower
     return,possible
end

