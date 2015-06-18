pro plot_res, fitsfile = fitsfile, datafile = datafile, log = log
device, decompose=0
if(NOT(keyword_set(datafile))) then datafile = '../output.data'
if(NOT(keyword_set(fitsfile))) then fitsfile = '../output.fits'

tab=readtab(datafile)

; first row
nuv= fix(tab[0, 0])
nv2= fix(tab[1, 0])
nt3amp = fix(tab[2, 0])
nvisamp = fix(tab[3, 0])
; second row
nt3phi = fix(tab[0, 1])
nt3 = fix(tab[1, 1])
nvisphi = fix(tab[2, 1])
nchanr =  fix(tab[3, 1])
; u,v,lambda
u=       reform(tab[0, 2:2+nuv-1])
v=       reform(tab[1, 2:2+nuv-1])
lambda = reform(tab[2, 2:2+nuv-1])
; t3 to uv index
t3in1 = reform(fix(tab[0, 2+nuv:2+nuv+nt3-1]))
t3in2 = reform(fix(tab[1, 2+nuv:2+nuv+nt3-1]))
t3in3 = reform(fix(tab[2, 2+nuv:2+nuv+nt3-1]))
; data & obs & residuals
reconst = reform(tab[0,2+nuv+nt3:* ])
data=    reform(tab[1, 2+nuv+nt3:*])
data_err= 1./reform(tab[2, 2+nuv+nt3:*])
res     = reform(tab[3, 2+nuv+nt3:*])

baseline = sqrt(u*u+v*v)

plotsym, 0, 0.5, /fill
loadct, 39
white = 0
blue = 64
yellow = 196
green = 128
orange = 220
red = 250

window, 0, xs=1000, ys=1000
!p.multi = [0, 1, 2] 
v2x = baseline[0:nv2-1]
v2y = data[0:nv2-1]
v2e = data_err[0:nv2-1]
v2y2 = reconst[0:nv2-1]
v2r = res[0:nv2-1]
if(NOT(keyword_set(log))) then plot, v2x, v2y, psym=8, title='V2 -- reconst & data vs baseline', xtitle='Spatial frequency (mega-lambda)', /xs, /ys else plot, v2x, v2y, psym=8, title='V2 -- reconst & data vs baseline', xtitle='Spatial frequency (mega-lambda)', /xs, /ys, /ylog, yr=[1e-4, 1]
oplot, v2x, v2y, col = blue, psym=8
errplot, v2x, v2y-v2e*0.5, v2y+v2e*0.5, col = blue
oplot, v2x, v2y2, col=red, psym=8
legend, ['Reconst', 'Data'], li=[0, 0], col=[red, blue], /right
plot, v2x, v2r,  psym=8, title='V2 -- Residuals', xtitle='Spatial frequency (mega-lambda)', ytitle='Normalized residuals', /xs


; TODO: check whether to use nt3 or nt3phi
window, 1, xs=1000, ys=1000
t3baseline = dblarr(nt3)
for i=0, nt3-1 do t3baseline[i] = min([baseline[t3in1[i]], baseline[t3in2[i]], baseline[t3in3[i]]])
t3phi_offset = nv2+nt3amp+nvisamp
!p.multi = [0, 1, 2] 
t3phix = t3baseline[0:nt3phi-1]
t3phiy = data[t3phi_offset:t3phi_offset+nt3phi-1] / !pi * 180.
t3phie = data_err[t3phi_offset:t3phi_offset+nt3phi-1]  / !pi * 180.
t3phiy2 = reconst[t3phi_offset:t3phi_offset+nt3phi-1] / !pi * 180.
t3phir = res[t3phi_offset:t3phi_offset+nt3phi-1] 
plot, t3phix, t3phiy, psym=8, title='Closure phases -- reconst & data vs baseline', xtitle='Spatial frequency (mega-lambda)', /xs, /ys, /nodata
oplot, t3phix, t3phiy, col = blue, psym=8
errplot, t3phix, t3phiy-t3phie*0.5, t3phiy+t3phie*0.5, col = blue, psym=8
oplot, t3phix, t3phiy2, col=red, psym=8
legend, ['Reconst', 'Data'], li=[0, 0], col=[red, blue], /right
plot, t3phix, t3phir,  psym=8, title='Closure phases -- Residuals', xtitle='Spatial frequency (mega-lambda)', ytitle='Normalized residuals', /xs

window, 2, xs=1000, ys=1000
!p.multi = [0, 1, 2] 
plot, v2y, v2y2, psym = 8, xtitle = 'Data', ytitle = 'Reconstructed', title='V2 -- Reconst vs data', /xs, /ys
plot, t3phiy, t3phiy2, psym = 8, xtitle = 'Data', ytitle = 'Reconstructed', title='Closures -- Reconst vs data', /xs, /ys


window, 3
loadct, 3
tab_im =  readfits(fitsfile,  head, /silent)
sz = size(tab_im)
if(sz[0] LT 2) then print, 'Error: image size is weird' else begin
if(sz[0] EQ 2) then nchanr = 1 else nchanr = (size(tab_im))[3]
ndf = 1.0*sxpar(head,  'NDF')
newchi2 = sxpar(head,  'CHISQR')
scale = sxpar(head,'scale')
xv=findgen(sz[1])*scale
yv=findgen(sz[2])*scale
xv=xv-mean(xv) 
xv=-xv
yv=yv-mean(yv)

if(nchanr LE 10) then begin
 window, 3, xs = 500*nchanr, ys = 500, title = 'Reconstructed Mean image'
       wset, 3
       multi_im=[0,nchanr,0] 
    endif else begin
       window, 3, xs = 500, ys = 500, title = 'Reconstructed Mean image'
       wset, 3
       multi_im=[0,1,0]
    endelse

!p.multi=multi_im
for i = 0, nchanr-1 do image_cont_uv, (reform(tab_im[*,*,i])),xv=xv,yv=yv,tit="Reconstruction - Channel"+strcompress(string(i)),xtit="RA (mas)",ytit="DEC (mas)", /asp,/noc

endelse

stop
end




plot_res
end
