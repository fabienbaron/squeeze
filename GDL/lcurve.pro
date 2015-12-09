pro lcurve ; simple lcurve example


; simulation part
file = '2004-data1'
tv_weight = [10, 100, 200, 500, 1000, 2000, 3000, 5000]
nfiles = n_elements(tv_weight)

for i=0, nfiles-1 do begin
command = '../bin/squeeze ../sample_data/'+file+'.oifits -tv '+strcompress(string(tv_weight[i]),/remove_all)+' -s 0.2 -w 70 -n 250 -e 2000 -o ../'+file+'_regtv'+strcompress(string(i),/remove_all)+'.fits'
spawn, command
endfor



; analysis part
chi2 = dblarr(nfiles) ; total reduced chi2 -- tied to neg-loglikelihood
tv_val = dblarr(nfiles) ; regularizer value

for i=0, nfiles-1 do begin
   im_name = '../2004-data1_regtv'+strcompress(string(i), /remove_all)+'.fits'
   im = readfits(im_name, head)
   tv_weight[i] =  sxpar(head, 'HYPER5')
   tv_val[i] =  sxpar(head, 'REGUL5W0')
   chi2[i] =  sxpar(head, 'CHI2')
;   image_cont, im, /asp, /noc, tit='a='+strcompress(string(tv_reg[i], format='(F0.1)'),/remove_all)
endfor

plot, tv_weight*tv_val, chi2, psym =-1, xtitle = 'Jprior', ytitle='Jdata', /xlog, /ylog, /xs, /ys
stop
end

lcurve
end
