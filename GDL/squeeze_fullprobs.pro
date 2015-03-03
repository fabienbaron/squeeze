; note: temperatures are only useful for parallel tempering
;
;
pro squeeze_fullprobs, file
device, decomposed=0
plotsym, 0, /fill
file = read_csv(file)

nthreads = fix(file.field1[0])
niter =  fix(file.field2[0])

; the number of temperatures is nthreads
ntemp_rows = nthreads / 3 + ((nthreads MOD 3) NE 0)
for i=1, ntemp_rows do begin
if(i EQ 1) then temp = [file.field1[i], file.field2[i], file.field3[i]] else temp = [temp, [file.field1[i], file.field2[i], file.field3[i]]]
endfor
temp = temp[0:nthreads-1]

like = reform(file.field1[ntemp_rows+1:*], niter, nthreads)
prior = reform(file.field2[ntemp_rows+1:*], niter, nthreads)
post = reform(file.field3[ntemp_rows+1:*], niter, nthreads)
loadct, 14


like_avg = dblarr(nthreads)
sigma_avg = dblarr(nthreads)
contrib = dblarr(nthreads)  
contrib_e = dblarr(nthreads)
;dumb averaging
off = ceil(0.5*niter)
for i=0, nthreads-1 do like_avg[i] = mean(reform(like[off:*, i]))
for i=0, nthreads-1 do sigma_avg[i] = (stdev(reform(like[off:*, i])))^2.

like_avg  = like_avg[sort(temp)]
temp = temp[sort(temp)]

print,  'Thread #      ', 'Temperature       ', 'Avg negloglikelihood  ', 'Stddev'
for i=0, nthreads-1 do print, i, temp[i], like_avg[i], sqrt(sigma_avg[i])
window, 0
plot, temp, like_avg, /xs, /ys, psym=8, /xlog, /ylog, xtit='Temperature', ytit='neglogLikelihood'
errplot, temp, like_avg-0.5*sqrt(sigma_avg), like_avg+0.5*sqrt(sigma_avg)

window, 1, xs=1000, ys=1000
loadct, 12
plot, like[*,0], yr=[100, 5000], col = 255, xtit='Iterations', ytit='neglogLikelihood'
oplot, [off, off], [100, 5000], col = 255, li =2, thick=4
for i=1, nthreads-1 do oplot, like[*, i], col = 255-(255/nthreads)*i
legend, string(temp), col = 255-(255/nthreads)*indgen(nthreads), li=0

logZ = 0
logZ_err = 0

print, 'Thread #  ', 'Temperature scaling  ', 'Avg negloglike i/i+1   ', 'LogZ contrib i/i+1 ',  'LogZe^2 contrib i-1/i+1'
for i = 0, nthreads-2 do begin
   contrib[i] = (1./temp[i+1]-1./temp[i])* 0.5 *(like_avg[i] + like_avg[i+1])
   if(i EQ 0) then contrib_e[i] = (sigma_avg[0]*(1./temp[1]-1./temp[0])^2. + sigma_avg[nthreads-1]*(1./temp[nthreads-1]-1./temp[nthreads-2])^2.) else contrib_e[i]= sigma_avg[i]*(1./temp[i+1]-1./temp[i-1])^2.
   logZ += contrib[i]
   logZ_err += contrib_e[i]
   print, i, (1./temp[i]-1./temp[i+1]), 0.5*(like_avg[i] + like_avg[i+1]), 0.5*contrib[i], 0.5*contrib_e[i]
endfor
logZ_err = sqrt(0.5*logZ_err)
print, 'logZ = ', logZ, '+/-', logZ_err


;newcontrib = contrib*0.
;newlogZ=0
;for i = 0, nthreads-2 do begin
;   newcontrib[i] = alog(exp(1./temp[i+1]-1./temp[i])* like_avg[i])
;   if(i EQ 0) then newcontrib_e[i] = (sigma_avg[0]*(1./temp[1]-1./temp[0])^2. + sigma_avg[nthreads-1]*(1./temp[nthreads-1]-1./temp[nthreads-2])^2.) else contrib_e[i]= sigma_avg[i]*(1./temp[i+1]-1./temp[i-1])^2.
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
