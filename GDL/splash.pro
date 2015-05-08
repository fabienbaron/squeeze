function splash, refimage, currentimage, metric = metric, range = range, subrange = subrange, dx = dx, dy = dy
;
; IN (mandatory):  refimage, currentimage: filenames of FITS files to compare 
; IN (optional) :  metric  : string, algorithm used for registration, e.g. 'mse', 'ssim'
;                  range   : integer, number of pixels defining the half-range searched for registration
;                  subrange: float, fraction of pixel to which the registration is done
;
; OUT           :   metric : value of the metric when images are co-registered
; OUT (optional):  dx, dy  : final shift
;
;

if (keyword_set(refimage) EQ 0) or (keyword_set(currentimage) EQ 0) then begin
print, 'ERROR in splash.pro : Missing image file' 
stop
return, 0
end

if keyword_set(metric) EQ 0 THEN metric = 'mse' ELSE BEGIN
indx = where(metric EQ ['mse', 'psnr', 'ssim', 'msssim', 'iqi'])
if(indx[0] EQ -1) then metric = 'mse'
ENDELSE

if(keyword_set(range) EQ 0) then range = 0
if(keyword_set(subrange) EQ 0) then subrange = 0

command = '/home/baron/SOFTWARE/splash/bin/splash --image1 '+refimage+' --image2 '+currentimage+' --algorithm '+metric+' --range '+string(range) +' --subrange '+string(subrange)+ ' > splash.tmp'
spawn, command
ret = (read_ascii('splash.tmp')).field1
dx = ret[4, 1]
dy = ret[5,1 ]
return, ret[3,1]
end
