pro idlmovie2

; ------ input -------------------------------------------------

 dimx0 = 400
 dimv0 = 400
 xmax0 = 5.0
 vmax0 = 1.0

 cutx  = 0
 dimx  = 400
 cutv  = 0
 dimv  = dimv0 - 2*cutv

 file01     = "phasex-1-sp0-"
 file02     = "phasex-2-sp0-"

 file_begin = 0.0
 file_end   = 7.0
 increment  = 1.0

; --------------------------------------------------------------

 special_file = -100

 s1     = dimx
 s2     = dimv
 xmin   = xmax0/dimx0 * cutx
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

		file1 = file01 + s 
		print, file1
		openu,1,file1
		a1 = assoc(1,bytarr(dimx0,dimv0,/nozero))
		phase1 = extrac(a1(0),cutx,cutv,dimx,dimv)

		file2 = file02 + s 
		print, file2
		openu,2,file2
		a2 = assoc(2,bytarr(dimx0,dimv0,/nozero))
		phase2 = extrac(a2(0),cutx,cutv,dimx,dimv)

		tvscl,alog10(1+phase1+phase2),xoff,yoff
		close,1
		close,2
	endrep until ( f gt file_end-increment )
 endif else begin
        if (special_file lt 10) then s = string( format='(F5.3)',special_file )
	if (special_file lt 0) then s = string( format='(F6.3)',special_file )
        if (special_file gt 10) then s = string( format='(F6.3)',special_file )
	file = file1 + s 
	print, file
	openu,1,file
	a = assoc(1,bytarr(dimx0,dimv0,/nozero))
	phase = extrac(a(0),cutx,cutv,dimx,dimv)
	tvscl,alog10(1+phase),xoff,yoff
	close,1
 endelse
end









