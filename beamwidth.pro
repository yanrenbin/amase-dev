function beamwidth,w1,npixel=npixel, fiber=fiber, fcoll=fcoll,sampling=sampling,alpha=alpha,d1=d1,refindex=refindex,halftheta=halftheta
; d1 is the distance between the prism outer surface and the entrace pupil of the camera in mm.

    if NOT keyword_set(fiber) then fiber=50.
    if NOT keyword_set(fcoll) then fcoll=105.
    if NOT keyword_set(sampling) then sampling=2.5
    if NOT keyword_set(npixel) then npixel=4096.
    if NOT keyword_set(refindex) then refindex=1.5

    w1arr = findgen(91)+10.
    halftheta = atan(npixel/(2.*sampling)*fiber/1000./fcoll) ; half angle outside the prism
    theta1 = asin(sin(halftheta)/refindex)  ; half of the angle within the prism
    w2arr = (w1arr*tan(alpha)*tan(theta1)+d1*tan(halftheta)+w1arr/2.)*2.
;    print,halftheta*!radeg
    plot,w1arr,w2arr/w1arr
    print,'Half angle outside the prism:',halftheta*!radeg
    return,interpol(w2arr,w1arr,w1)
end
