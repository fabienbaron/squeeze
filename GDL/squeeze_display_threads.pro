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

; This is a GDL viewer for SQUEEZE that constantly monitor the various threads
;
;  keyword TEMPSORT: forces sorting by temperatures, useful for parallel tempering

pro squeeze_threads, dir, tempsort = tempsort

dir = dir +'/'
if not(keyword_set(dir)) then dir ='./'

device, decompose=0
loadct, 3
restart:
ltime =  systime(1)
mtime =  ltime
while (mtime-ltime le  0) do begin
   openr,  1,  'thread00.fits',  error = err
   if (err eq 0) then begin
      stat =  fstat(1)
      close,  1
      mtime =  stat.mtime
   endif else print,  'Trouble opening thread00.fits (ignore if occasional message only)'
   wait,  0.3
endwhile
ltime = mtime
tab_im =  readfits('thread00.fits',  head, /silent)
sz = size(tab_im)
if(sz[0] LT 2) then goto, restart
nthreads = sxpar(head,  'NTHREADS')


!p.multi=[0,nthreads,1]
window, 0, xs = 250*nthreads, ys = 250, title = 'Threads', xpos = 0, ypos = 0


temperatures = dblarr(nthreads)
images = PTRARR(nthreads, /allocate_heap)
while (1) do begin

   for i=0, nthreads-1 do begin
      filename = 'thread'+strcompress(string(i, FORMAT='(I02)'),/remove_all)+'.fits'
      if( FILE_TEST(filename, /read) EQ 1) then begin
         *images[i] = readfits(filename, head, /silent)
         temperatures[i] = sxpar(head,'TEMPER') 
      endif
   endfor 
   reorder = sort(temperatures)

   for i=0, nthreads-1 do begin
         if(keyword_set(tempsort)) then j = reorder[i] else j = i
         im = *images[j]
         sz = size(im)
         if(sz[0] GE 2) then image_cont_uv, sqrt(reform(im)),tit="Thread "+strcompress(string(j))+" - temp "+strcompress(string(temperatures[j])),xtit="pixels",ytit="pixels", /asp,/noc else image_cont_uv, dblarr(32, 32), tit="Thread "+strcompress(string(j))+" - temp "+strcompress(string(temperatures[j])),xtit="pixels",ytit="pixels", /asp,/noc     
      endfor
wait, 3
endwhile

end


; TESTING EXAMPLES
;./squeeze ./sample_data/2004-data1.fits -w 128 -s 0.2 -threads 10 -i random -n 2000
