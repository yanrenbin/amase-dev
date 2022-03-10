pro bkgd_brightness

   list=file_search('~/astro/stellib/platedesign/2019future/gaia/','*_all-result.fits.gz',count=ct0)
   k=where((strmatch(list,'*+??_btx*.fits.gz') or strmatch(list,'*-??_btx*.fits.gz')) and strmatch(list,'*SGR1-RV_btx*.fits.gz') eq 0,ct)

   pos=strpos(list[k],'_btx')
   lat = intarr(ct)
   for i=0,ct-1 do lat[i]=fix(strmid(list[k[i]],pos[i]-2,2))
   t = where(lat eq 0)
   lat[t]=4
   nanoperarcsecsq = fltarr(ct)
   for i=0,ct-1 do begin
      obj=mrdfits(list[k[i]],1)
      nanomaggies=10^(-0.4*(obj.phot_g_mean_mag-22.5))
      nanoperarcsecsq[i] = total(nanomaggies)/(!pi*1.49^2)/3600./3600.
   endfor
   openw,lun,'bkgd_brightness.txt',/get_lun
   for i=0,ct-1 do printf,lun,lat[i],nanoperarcsecsq[i],format='(i,f)'
   close,lun
   free_lun,lun
end

