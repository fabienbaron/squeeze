;
; SQUEEZE utility, reads and display the full MCMC chain
;

pro squeeze_fullchain, filename

; dimensions are axis_len, axis_len, nchanr, niter, nchains
device, decomposed = 0
loadct, 3
images = readfits(filename, head)
sz = size(images)
nchains = sxpar(head,'nchains') 
niter = sxpar(head,'niter')  
nchanr   = sxpar(head,'nchanr')  
nelements = sxpar(head,'elements')  
axis_len = sz[0]  

print, 'Nchains: ', nchains 
print, 'Niter: ', niter
print, 'Nchanr: ', nchanr
print, 'Nelements: ', nelements
print, 'Image width: ', axis_len

snplots = ceil(sqrt(niter))

for t=0, nchains-1 do begin
   window, t, xs = 50 * snplots, ys = 50 * snplots, title = 'Chain'+string(t)
   !p.multi = [0, snplots, snplots]
   !x.margin = 0
   !y.margin = 0
   !z.margin = 0
   for i = 0, niter-1 do image_cont_uv, sqrt(reform(images[*,*, 0, i, t])), /asp, /noc
endfor

stop

end
