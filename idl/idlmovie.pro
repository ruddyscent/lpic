pro idlmovie

 file0      = "phasex-sp0-"
 file_begin = 0.0
 file_end   = 10.0
 increment  = 1.0

 special_file = -100

 dimx0 = 400
 dimv0 = 400
 xmax0 = 5.0
 vmax0 = 1.0

 cutx  = 0
 dimx  = 400
 cutv  = 0
 dimv  = dimv0 - 2*cutv

; ----------------------------------------------------------------------------------------

 s1     = dimx
 s2     = dimv
 xmin   = xmax0/dimx0 * cutx - xoffset
 xmax   = xmin + xmax0/dimx0 * dimx
 ymax   = vmax0 - vmax0/dimv0 * cutv  
 x      = xmin + (xmax-xmin)/s1 * findgen(s1)
 y      = -ymax + 2.0*ymax/s2 * findgen(s2)

 print, cutv
 print, ymax

 set_plot,"x"
 window, 1, xpos=540, ypos=28, xsize=dimx+200, ysize=dimv+180, title="phasespace(x,v)"
 xoff=120
 yoff=100

 loadct,13
 tvlct, r, g, b, /get
 r(0) = 255
 g(0) = 255
 b(0) = 255
 r(1) = 0
 g(1) = 0
 b(1) = 0
 tvlct, r, g, b
 	 
 plot, x, y, xstyle=1, ystyle=1, xrange=[xmin,xmax], yrange=[-ymax,ymax],$
	xtitle="!3x/!7k!3", ytitle="!3v/c!3",ticklen=-0.02,$
	charsize=1.5, color=1, position=[xoff-1,yoff-1,xoff+s1+1,yoff+s2+1],$
	/noerase, /nodata, /device

 if (special_file lt -10) then begin
	i=-1
	repeat	begin
		i = i+1
		f = floor( 1000.0*(file_begin+increment*i)+0.5 ) / 1000.0
		if (f lt 10) then s = string( format='(F5.3)',f )
		if (f ge 10) then s = string( format='(F6.3)',f )
		file = file0 + s 
		print, file
		openu,1,file
		a = assoc(1,bytarr(dimx0,dimv0,/nozero))
		phase = extrac(a(0)+1,cutx,cutv,dimx,dimv)
		tvscl,alog10(phase),xoff,yoff
		close,1
	endrep until ( f gt file_end-increment )
 endif else begin
        if (special_file lt 10) then s = string( format='(F5.3)',special_file )
	if (special_file lt 0) then s = string( format='(F6.3)',special_file )
        if (special_file gt 10) then s = string( format='(F6.3)',special_file )
	file = file0 + s 
	print, file
	openu,1,file
	a = assoc(1,bytarr(dimx0,dimv0,/nozero))
	phase = extrac(a(0)+1,cutx,cutv,dimx,dimv)
	tvscl,alog10(phase),xoff,yoff
	close,1
 endelse
end









