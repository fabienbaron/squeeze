FUNCTION readnextline,lun,line
  ok=0
  while ( (not eof(lun)) and (ok eq 0) ) do begin
     readf,lun,line
     bl= byte(strtrim(line,2))
     ok=1                     
     if  ( strlen(line) eq 0 ) or ( (bl(0) EQ 35) OR (bl(0) EQ 59) OR (bl(0) EQ 47) ) then begin
        ok=0                    ; lines beginning by ; or / or # are comments
     endif else begin
                                ; modify line if it contains text
        bl2 = bl * ( ( (bl GE 43) * (bl LE 57)  ) OR (bl EQ 32) OR (bl EQ 9) OR (bl EQ 101) OR (bl EQ 69))
;stop
        indx = where(bl2 EQ 0)
        if(indx[0] NE -1) then begin
           bl2[where(bl2 EQ 0)] = byte('0')
           line = string(bl2)
        endif
     endelse
  end
  return,ok
end


FUNCTION readtab, filename, cols=cols

  print, 'Reading as txt file ', filename
  openr,lun,filename,/get_lun
  line=''
  if readnextline(lun,line) eq 0 then message,'No data'
  bl=byte(line)
  
  s= size(bl) & n= s[1] - 1
  ncol=0L
  ws=0L
  wsold=1L
  ii=0L
  while (ii LE n) do begin
      if (bl[ii] eq 32 or bl[ii] eq 09 ) then ws=1 else ws=0 ; 
      if ws ne wsold then begin
          if (ws eq 0) then ncol=ncol+1
	  wsold=ws
      end
  ii=ii+1L
  end

  print,'File contains', ncol,' columns'

  datline= dblarr(ncol)
  point_lun,lun,0
  ii=0L
  while not eof(lun) do begin
      if readnextline(lun,line) eq 1 then begin
          ii=ii+1
          reads,line,datline
          if ii eq 1 then data=datline  else data= [[data],[datline]]
      end
  end
 
  close,lun & free_lun,lun

  if keyword_set(cols) then begin
      cols=cols-1   
      data=data([cols],*)
  end
  return,data
end






