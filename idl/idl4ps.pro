pro idl4ps

 file1  = "pic.phasex_sp2.25.000"
 title1 = "phasespace(x,vx) t=25.0"

 file2  = "pic.phasex_sp2.51.000"
 title2 = "phasespace(x,vx) t=51.0"

 file3  = "pic.phasex_sp2.75.000"
 title3 = "phasespace(x,vx) t=75.0"

 file4  = "pic.phasex_sp2.99.000"
 title4 = "phasespace(x,vx) t=99.0"

 dimx0 = 400
 dimv0 = 400
 xmax0 = 20.46
 vmax0 = 0.5

 cutx  = 0
 dimx  = 400
 cutv  = 0
 dimv  = 400 - 2*cutv

 xmax   = 20.46
 ymax   = 0.5

 x      = xmax/dimx * findgen(dimx)
 y      = -ymax + 2.0*ymax/dimv * findgen(dimv)
 xcm    = 29
 ycm    = 21

 set_plot,"ps"

 xpix = xcm * !d.x_px_cm
 ypix = ycm * !d.y_px_cm
 xoff = !d.x_px_cm * 2.5 
 yoff = !d.x_px_cm * 4.0
 sx   = !d.x_px_cm * 7.5 
 sy   = !d.y_px_cm * 8.0
 dx   = xoff + sx 
 dy   = !d.y_px_cm * 11.0

 device, file="phase.ps", /color, bits=8, /portrait,$
         xoffset=0, yoffset=0, xsize=xpix, ysize=ypix
  
 xo = xoff 
 yo = yoff + dy
 openu,1,file1
 a = assoc(1,bytarr(dimx0,dimv0,/nozero))
 phase = extrac(a(0)+1,cutx,cutv,dimx,dimv)
 tvscl, alog10(phase), xsize=sx, ysize=sy, xo, yo
 plot, x, y, xstyle=1, ystyle=1, xrange=[0,xmax], yrange=[-ymax,ymax],$
       xtitle="!3x/!7k!3", ytitle="!3vx/c!3", title=title1,$
       color=1, position=[xo,yo,xo+sx,yo+sy], /noerase, /nodata, /device 
 close,1 

 xo = xoff + dx
 yo = yoff + dy
 openu,1,file2
 a = assoc(1,bytarr(dimx0,dimv0,/nozero))
 phase = extrac(a(0)+1,cutx,cutv,dimx,dimv)
 tvscl, alog10(phase), xsize=sx, ysize=sy, xo, yo
 plot, x, y, xstyle=1, ystyle=1, xrange=[0,xmax], yrange=[-ymax,ymax],$
       xtitle="!3x/!7k!3", ytitle="!3vx/c!3", title=title2,$
       color=1, position=[xo,yo,xo+sx,yo+sy], /noerase, /nodata, /device 
 close,1

 xo = xoff 
 yo = yoff
 openu,1,file3
 a = assoc(1,bytarr(dimx0,dimv0,/nozero))
 phase = extrac(a(0)+1,cutx,cutv,dimx,dimv)
 tvscl, alog10(phase), xsize=sx, ysize=sy, xo, yo
 plot, x, y, xstyle=1, ystyle=1, xrange=[0,xmax], yrange=[-ymax,ymax],$
       xtitle="!3x/!7k!3", ytitle="!3vx/c!3", title=title3,$
       color=1, position=[xo,yo,xo+sx,yo+sy], /noerase, /nodata, /device 
 close,1

 xo = xoff + dx
 yo = yoff 
 openu,1,file4
 a = assoc(1,bytarr(dimx0,dimv0,/nozero))
 phase = extrac(a(0)+1,cutx,cutv,dimx,dimv)
 tvscl, alog10(phase), xsize=sx, ysize=sy, xo, yo
 plot, x, y, xstyle=1, ystyle=1, xrange=[0,xmax], yrange=[-ymax,ymax],$
       xtitle="!3x/!7k!3", ytitle="!3vx/c!3", title=title4,$
       color=1, position=[xo,yo,xo+sx,yo+sy], /noerase, /nodata, /device 
 close,1

 device,/close
 
 set_plot, "x"

end



