; $Id: image_cont.pro,v 1.1 1993/04/02 19:43:31 idl Exp $

pro image_cont_uv, a, WINDOW_SCALE = window_scale, ASPECT = aspect, $
	INTERP = interp, nocontours=nocontours,xval=x,yval=y,xtit=xtit,$
        ytit=ytit,noimage=noimage,levs=levs,cann=cann,$
        tit=tit,subtit=subtit,possubtit=possubtit,logim=logim,colsty=colsty

;+
; NAME:
;	IMAGE_CONT
;
; PURPOSE:
;	Overlay an image and a contour plot.
;
; CATEGORY:
;	General graphics.
;
; CALLING SEQUENCE:
;	IMAGE_CONT, A
;
; INPUTS:
;	A:	The two-dimensional array to display.
;
; KEYWORD PARAMETERS:
; WINDOW_SCALE:	Set this keyword to scale the window size to the image size.
;		Otherwise, the image size is scaled to the window size.
;		This keyword is ignored when outputting to devices with 
;		scalable pixels (e.g., PostScript).
;
;	ASPECT:	Set this keyword to retain the image's aspect ratio.
;		Square pixels are assumed.  If WINDOW_SCALE is set, the 
;		aspect ratio is automatically retained.
;
;	INTERP:	If this keyword is set, bilinear interpolation is used if 
;		the image is resized.
;       LOGIM : pass the cut off value and this will
;		overplot the logarithm of the data above cutoff.
;      COLSTY : Color Style - set to different values depending upon which 
;               lookup table to use for color
;
; OUTPUTS:
;	No explicit outputs.
;
; COMMON BLOCKS:
;	None.
;
; SIDE EFFECTS:
;	The currently selected display is affected.
;
; RESTRICTIONS:
;	None.
;
; PROCEDURE:
;	If the device has scalable pixels, then the image is written over
;	the plot window.
;
; MODIFICATION HISTORY:
;	DMS, May, 1988.
;       JDM, 1996.       Added nocontour keyword, and 
;  			 changed tv--> tvscl in one subroutine
; 	JDM, 1997.       (+_deluxe) Snazzed it up for 
;			 Publication Quality Use.
;       PGT, 1998 Nov    Added colsty keyword for better control of colors
;	JDM, 2003Jul	 should allow inversion of axes if xv,yv are inverted.
;			 will not invert image, only label is right!
;-

on_error,2                      ;Return to caller if an error occurs
sz = size(a)			;Size of image
if (keyword_set(x) eq 0) then x=findgen(sz(1))
if (keyword_set(y) eq 0) then y=findgen(sz(2))
if (keyword_set(tit) eq 0) then tit=' '
if (keyword_set(xtit) eq 0) then xtit=' '
if (keyword_set(ytit) eq 0) then ytit=' '
if (keyword_set(noimage) eq 0) then noimage=0
if (keyword_set(logim) eq 0) then logim=0.0

if (keyword_set(levs) eq 0) then levs=0
if (keyword_set(cann) eq 0) then cann=0
if (keyword_set(colsty) eq 0) then colsty=0
 
if sz(0) lt 2 then message, 'Parameter not 2D'

	;set window used by contour
contour,[[0,0],[1,1]],/nodata, xstyle=4, ystyle = 4

px = !x.window * !d.x_vsize	;Get size of window in device units
py = !y.window * !d.y_vsize
swx = px(1)-px(0)		;Size in x in device units
swy = py(1)-py(0)		;Size in Y
six = float(sz(1))		;Image sizes
siy = float(sz(2))
aspi = six / siy		;Image aspect ratio
aspw = swx / swy		;Window aspect ratio
f = aspi / aspw			;Ratio of aspect ratios

if (!d.flags and 1) ne 0 then begin	;Scalable pixels?
  if keyword_set(aspect) then begin	;Retain aspect ratio?
				;Adjust window size
	if f ge 1.0 then swy = swy / f else swx = swx * f
	endif
if (keyword_set(noimage) eq 0) then begin
    if (logim eq 0.0) then begin
        if (colsty eq 0) then begin
          tvscl,a,px(0),py(0),xsize = swx, ysize = swy, /device
        endif
        if (colsty eq 1) then begin
           tmp=(a-min(a))/(max(a)-min(a))*194. + 1
           tv,tmp,px(0),py(0),xsize = swx, ysize = swy, /device
        endif
        if (colsty eq 2) then begin
           tmp=(a-min(a))/(max(a)-min(a))*195.
           tv,tmp,px(0),py(0),xsize = swx, ysize = swy, /device
        endif
    endif else begin
        if (colsty eq 0) then begin
          tvscl,-1*alog(a>logim),px(0),py(0),xsize = swx, ysize = swy, /device
        endif 
        if (colsty eq 1) then begin
           tmp=alog(a>logim)
           tmp=(tmp-min(tmp))/(max(tmp)-min(tmp))*194. + 1
           tv,tmp,px(0),py(0),xsize = swx, ysize = swy, /device
        endif 
        if (colsty eq 2) then begin
           tmp=alog(a>logim)
           tmp=(tmp-min(tmp))/(max(tmp)-min(tmp))*195.
           tv,tmp,px(0),py(0),xsize = swx, ysize = swy, /device
        endif
    endelse
   endif


  

endif else begin	;Not scalable pixels	
   if keyword_set(window_scale) then begin ;Scale window to image?

if (keyword_set(noimage) eq 0) then begin
    if (logim eq 0.0) then begin
	tvscl,a,px(0),py(0)	;Output image
    endif else begin
        print,'location1'
        tvscl,alog(a>logim),px(0),py(0)
    endelse

      
endif
	swx = six		;Set window size from image
	swy = siy
    endif else begin		;Scale window
	if keyword_set(aspect) then begin
		if f ge 1.0 then swy = swy / f else swx = swx * f
		endif		;aspect
; I changed the following line to read 'tvscl,...' instead of 'tv,...'
;     --JDM 11/19/96
if (keyword_set(noimage) eq 0) then begin
  if (logim eq 0.0) then begin
	tvscl,poly_2d(bytscl(a),$	;Have to resample image
		[[0,0],[six/swx,0]], [[0,siy/swy],[0,0]],$
		keyword_set(interp),swx,swy), $
		px(0),py(0)
      endif
   if (logim ne 0.0) then begin 
         print,'location2'
         tvscl,poly_2d(bytscl(alog(a>logim)),$       ;Have to resample image
                [[0,0],[six/swx,0]], [[0,siy/swy],[0,0]],$
                keyword_set(interp),swx,swy), $
                px(0),py(0)
      endif
      endif
	endelse			;window_scale
  endelse			;scalable pixels

mx = !d.n_colors-1		;Brightest color
colors = [mx,mx,mx,0,0,0]	;color vectors
colors= [mx,mx,mx,mx,mx,mx]

if !d.name eq 'PS' then begin
  colors = mx - colors ;invert line colors for pstscrp
  
endif
if (keyword_set(nocontours) eq 0) then begin
in_p=where(levs ge 0,p_ct)
in_n=where(levs lt 0,n_ct)
if (p_ct ne 0) then begin
  plevs=levs(in_p)
  plevs=plevs(sort(plevs))
endif
if (n_ct ne 0) then begin
  nlevs=levs(in_n)
  nlevs=nlevs(sort(nlevs))
endif


; First do positive contours
if (keyword_set(subtit) eq 1) then begin
;Do this in Normalized Coords.
if (keyword_set(possubtit) eq 0) then possubtit=[.13,.05]
xyouts,possubtit(0),possubtit(1),charsize=.5,subtit,/norm
endif
if (p_ct gt 0) then $
contour,a,x,y,/noerase,/xst,/yst,$	;Do the contour
	   pos = [px(0),py(0), px(0)+swx,py(0)+swy],/dev,$
	c_color =  colors ,xtit=xtit,ytit=ytit,tit=tit,levels=plevs ,$
	/follow,c_lab=0,$
   xr=[x(0),x(n_elements(x)-1)], yr =[y(0),y(n_elements(y)-1)]
if (n_ct gt 0) then $
contour,a,x,y,/noerase,/xst,/yst,$      ;Do the contour
           pos = [px(0),py(0), px(0)+swx,py(0)+swy],/dev,$
        c_color =  colors ,xtit=xtit,ytit=ytit,tit=tit,levels=nlevs ,$
	c_line=2,/follow,c_lab=0,$
   xr=[x(0),x(n_elements(x)-1)], yr =[y(0),y(n_elements(y)-1)]

;print,'neg'



endif $
else contour,a,x,y,/noerase,/xst,/yst,$  ;Do the contour
           pos = [px(0),py(0), px(0)+swx,py(0)+swy],/dev,$
        c_color =  colors ,/nodata,xtit=xtit,ytit=ytit,tit=tit,levels=levs,$
   xr=[x(0),x(n_elements(x)-1)], yr =[y(0),y(n_elements(y)-1)]



return
end
