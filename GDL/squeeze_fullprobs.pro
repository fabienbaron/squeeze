; note: temperatures are only useful for parallel tempering
;
;
pro squeeze_fullprobs, file
device, decomposed=0
plotsym, 0, /fill
file = read_csv(file)

nchains = fix(file.field1[0])
niter =  fix(file.field2[0])

; the number of temperatures is nchains
ntemp_rows = nchains / 3 + ((nchains MOD 3) NE 0)
for i=1, ntemp_rows do begin
if(i EQ 1) then temp = [file.field1[i], file.field2[i], file.field3[i]] else temp = [temp, [file.field1[i], file.field2[i], file.field3[i]]]
endfor
temp = temp[0:nchains-1]

like = reform(file.field1[ntemp_rows+1:*], niter, nchains)
prior = reform(file.field2[ntemp_rows+1:*], niter, nchains)
post = reform(file.field3[ntemp_rows+1:*], niter, nchains)
loadct, 14


like_avg = dblarr(nchains)
sigma_avg = dblarr(nchains)
contrib = dblarr(nchains)  
contrib_e = dblarr(nchains)
;dumb averaging
off = ceil(0.5*niter)
for i=0, nchains-1 do like_avg[i] = mean(reform(like[off:*, i]))
for i=0, nchains-1 do sigma_avg[i] = (stdev(reform(like[off:*, i])))^2.

like_avg  = like_avg[sort(temp)]
temp = temp[sort(temp)]

print,  'Chain #      ', 'Temperature       ', 'Avg negloglikelihood  ', 'Stddev'
for i=0, nchains-1 do print, i, temp[i], like_avg[i], sqrt(sigma_avg[i])
window, 0
plot, temp, like_avg, /xs, /ys, psym=8, /xlog, /ylog, xtit='Temperature', ytit='neglogLikelihood'
errplot, temp, like_avg-0.5*sqrt(sigma_avg), like_avg+0.5*sqrt(sigma_avg)

window, 1, xs=1000, ys=1000
loadct, 12
plot, like[*,0], yr=[100, 5000], col = 255, xtit='Iterations', ytit='neglogLikelihood'
oplot, [off, off], [100, 5000], col = 255, li =2, thick=4
for i=1, nchains-1 do oplot, like[*, i], col = 255-(255/nchains)*i
legend, string(temp), col = 255-(255/nchains)*indgen(nchains), li=0

logZ = 0
logZ_err = 0

print, 'Chain #  ', 'Temperature scaling  ', 'Avg negloglike i/i+1   ', 'LogZ contrib i/i+1 ',  'LogZe^2 contrib i-1/i+1'
for i = 0, nchains-2 do begin
   contrib[i] = (1./temp[i+1]-1./temp[i])* 0.5 *(like_avg[i] + like_avg[i+1])
   if(i EQ 0) then contrib_e[i] = (sigma_avg[0]*(1./temp[1]-1./temp[0])^2. + sigma_avg[nchains-1]*(1./temp[nchains-1]-1./temp[nchains-2])^2.) else contrib_e[i]= sigma_avg[i]*(1./temp[i+1]-1./temp[i-1])^2.
   logZ += contrib[i]
   logZ_err += contrib_e[i]
   print, i, (1./temp[i]-1./temp[i+1]), 0.5*(like_avg[i] + like_avg[i+1]), 0.5*contrib[i], 0.5*contrib_e[i]
endfor
logZ_err = sqrt(0.5*logZ_err)
print, 'logZ = ', logZ, '+/-', logZ_err


;newcontrib = contrib*0.
;newlogZ=0
;for i = 0, nchains-2 do begin
;   newcontrib[i] = alog(exp(1./temp[i+1]-1./temp[i])* like_avg[i])
;   if(i EQ 0) then newcontrib_e[i] = (sigma_avg[0]*(1./temp[1]-1./temp[0])^2. + sigma_avg[nchains-1]*(1./temp[nchains-1]-1./temp[nchains-2])^2.) else contrib_e[i]= sigma_avg[i]*(1./temp[i+1]-1./temp[i-1])^2.
   newlogZ += contrib[i]
;   newlogZ_err += contrib_e[i]
;   print, i, (1./temp[i]-1./temp[i+1]), 0.5*(like_avg[i] + like_avg[i+1]), 0.5*contrib[i], 0.5*contrib_e[i]
;endfor
;print, 'other logZ expression: ', newlogZ



window,2 
plot, temp, contrib, psym = 8, tit='Log Z contribution', xtit='Temperature', ytit = 'Log Z contrib', /xlog

window, 3
plot, temp, contrib_e, psym = 8, tit='Log Ze contribution', xtit='Temperature', ytit = 'Log Ze contrib', /xlog

window, 4
plot, temp, like_avg, psym = 8, tit='Avg neglogLikelihood', xtit='Temperature', ytit = 'Avg neglogLikelihood', /xlog, /ylog
stop
end
