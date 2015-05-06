using FITSIO
using Tk
using ImageView
using Images
fullchain_filename = "../output_fullchain.fits.gz"
f = FITS(fullchain_filename)
data=read(f[1])

c = canvasgrid(16,16)
ops = [:pixelspacing => [1,1]]

for i=0:250
 img = data[:,:,1, i+1]
 ix = div(i,16)+1
 iy = i%16+1
 view(c[ix,iy], img)
end
readline()