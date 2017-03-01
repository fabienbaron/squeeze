pro recenter_imstack, im
; NxNxNWxNC image stack, i.e. NC chains of images of (NxN) pixels x NW wavelengths
m=(size(im))[1]
writefits, 'ref.fits', reform(im[0, *,*])
for t=1, m-1 do begin
   writefits, 'try.fits', reform(im[t, *,*])
   print, splash('ref.fits', 'try.fits', dx=dx, dy=dy, range = 20 , subrange = .1)
   print, dx, dy
   im[t, *,*] = fftshift(reform(im[t, *,*]), dx, dy)
endfor

end
