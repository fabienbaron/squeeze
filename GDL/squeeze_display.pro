;
;    SQUEEZE is free software: you can redistribute it and/or modify
;    it under the terms of the GNU General Public License as published by
;    the Free Software Foundation, either version 3 of the License, or
;    (at your option) any later version.
;
;    SQUEEZE is distributed in the hope that it will be useful,
;    but WITHOUT ANY WARRANTY; without even the implied warranty of
;    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
;    GNU General Public License for more details.
;
;    You should have received a copy of the GNU General Public License
;    along with Squeeze.  If not, see <http://www.gnu.org/licenses/>.

; This is a GDL viewer for SQUEEZE based on John Monnier's code for MACIM
;
pro squeeze_display, dir
if not(keyword_set(dir)) then dir ='../' else dir = dir +'/'
print, 'Monitoring chain00.fits in '+dir
device,retain=2, decompose=0
nreguls = 9
reg_names = ['PARAM', 'CENT', 'IMPRIOR', 'ENT', 'DEN', 'TV', 'SPOT', 'LAP', 'L0', 'TRANSPEC']
loadct,3
ltime =  systime(1)
mtime =  ltime
chi2 = -1
restart:                        ;restart from here if severe file read failure for chain.fits
while (1) do begin

; monitor file changes
 while(1) do begin 
 openr, 1,  dir+'chain00.fits',  error = err
  if (err eq 0) then begin
   stat =  fstat(1)
   close,  1
   mtime =  stat.mtime
   ;print, mtime
  endif else print,  'Trouble opening file (ignore if occasional message only)'
  if (mtime GT ltime) then break else wait, 0.3
 endwhile
 ltime = mtime		
 tab_im =  readfits(dir+'chain00.fits',  head, /silent)
 sz = size(tab_im)

 if(sz[0] LT 2) then goto, restart
 if(sz[0] EQ 2) then nchanr = 1 else nchanr = (size(tab_im))[3]
 if(n_elements(oldnchanr) EQ 0) then oldnchanr = nchanr
 
 if ( (n_elements(multi_im) EQ 0) OR (oldnchanr NE nchanr) ) then begin
    oldnchanr = nchanr
    if(nchanr LE 10) then begin
       window, 0, xs = 250*nchanr, ys = 500, title = 'Current/Mean image', xpos = 505
       wset, 0
       multi_im=[0,nchanr,2] 
    endif else begin
       window, 0, xs = 250, ys = 500, title = 'Current/Mean image', xpos = 505
       wset, 0
       multi_im=[0,1,2]
    endelse
 endif
 ndf = 1.0*sxpar(head,  'NDF')
 newchi2 = sxpar(head,  'CHISQR')
 scale = sxpar(head,'scale')
 if n_elements(reg_params) EQ 0 then reg_params = dblarr(nreguls)
 if n_elements(reg_vals) EQ 0 then reg_vals = dblarr(nchanr, nreguls)
 for i=0, nreguls-1 do begin
    reg_params[i] = sxpar(head,  'HYPER'+strcompress(string(i),/remove_all))
    if(reg_params[i] GT 0) then for j=0, nchanr-1 do reg_vals[j, i] = sxpar(head,  'REG'+strcompress(string(i),/remove_all)+'W'+strcompress(string(j),/remove_all))
 endfor

reg_active = where(reg_params GT 0)
valid = 0
if (n_elements(reg_active) GT 1) then valid = 1 else if ((n_elements(reg_active) EQ 1) AND (reg_active NE -1) ) then valid = 1
if valid EQ 1 then begin
   newregs = dblarr(nchanr, n_elements(reg_active))
   for i=0, nchanr-1 do newregs[i, *] = reg_vals[i, reg_active]*reg_params[reg_active]/(2.*ndf)
   linesty = (indgen(9))[reg_active]
endif


xv=findgen(sz[1])*scale
yv=findgen(sz[2])*scale
xv=xv-mean(xv) 
xv=-xv
yv=yv-mean(yv)
wset, 0
!p.multi=multi_im
for i = 0, nchanr-1 do image_cont_uv, sqrt(reform(tab_im[*,*,i])),xv=xv,yv=yv,tit="Current - Channel"+strcompress(string(i)),xtit="mas",ytit="mas", /asp,/noc


if FILE_TEST(dir+'output.fits') NE 1 then goto, skip

tab_mim=readfits(dir+'output.fits', head2,/silent)
szm = size(tab_mim) 
if(sz[0] LT 2) then goto, skip
if(szm[0] EQ 2) then nchanrm = 1 else nchanrm = szm[3]
scale=sxpar(head2,'scale')  
xv=findgen(szm[1])*scale
yv=findgen(szm[2])*scale
xv=xv-mean(xv)
xv=-xv
yv=yv-mean(yv)
if(nchanr EQ  nchanrm) then for i = 0, nchanrm-1 do image_cont_uv, sqrt(reform(tab_mim[*,*,i])),xv=xv,yv=yv,tit="Last (avg) - Channel"+strcompress(string(i)),xtit="mas",ytit="mas", /asp,/noc

skip:                           ; goto here if severe read failure for output.fits

if (chi2[0] EQ -1) then begin
   chi2 = newchi2
   regs = newregs
   window, 2, xs = 500 , ys = 500, xpos = 0, title = 'Chi2/regularizers vs reconstruction trials'
endif else begin
   wset, 2
   chi2 =  [chi2, newchi2]
   if((size(regs))[2] EQ (size(newregs))[2]) then regs = [regs, newregs] else regs = newregs
   !p.multi=[0,0,0]
   plot, chi2, ytit="Chi2r & Reguls/(2*ndf)",xtit="Trial",tit="Chi2 & Regularizations", yminor = 20, yr=[min(regs[*,reg_active]), max(chi2)], /ylog, /xs
   if valid EQ 1 then begin
  ; for j= 0, nchanr-1 do begin
      ;stop
      for i=0, n_elements(reg_active)-1 do oplot, regs[*,i], li=linesty[i];, col=j*64   
  ; endfor
   legend, ['CHI2', reg_names[reg_active]], li=[0,linesty]
endif
   wait, 1 
endelse

endwhile
finish:
print,  'Enter pressed, exiting...'

end
