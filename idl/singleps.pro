pro single2ps

 dimx0 = 400
 dimv0 = 400
 xmax0 = 4.095
 vmax0 = 1.0

 cutx  = 150
 dimx  = 200
 cutv  = 50
 dimv  = dimv0 - 2*cutv

 file      = "pic.phasex_sp1.10.200"

 s1     = dimx
 s2     = dimv
 xmin   = xmax0/dimx0 * cutx
 xmax   = xmin + xmax0/dimx0 * dimx
 ymax   = vmax0 - vmax0/dimv0 * cutv  
 x      = xmin + (xmax-xmin)/s1 * findgen(s1)
 y      = -ymax + 2.0*ymax/s2 * findgen(s2)

 set_plot,"x"
 window, 1, xpos=540, ypos=28, xsize=dimx+200, ysize=dimv+200, title="phasespace(x,v)"
 xoff=100
 yoff=100

 print, file
 openu,1,file
 a = assoc(1,bytarr(dimx0,dimv0,/nozero))
 phase = extrac(a(0)+1,cutx,cutv,dimx,dimv)
 tvscl,alog10(phase),xoff,yoff
 close,1

 plot, x, y, xstyle=1, ystyle=1, xrange=[xmin,xmax], yrange=[-ymax,ymax],$
	xtitle="!3x/!7k!3", ytitle="!3vx/c!3",ticklen=0.02,$
	color=1, position=[xoff-1,yoff-1,xoff+s1+1,yoff+s2+1],$
	/noerase, /nodata, /device

end









