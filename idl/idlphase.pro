pro idlphase

; ------ input ---------------------------------------------------------------------------

 dimx0 = 400
 dimv0 = 400
 xmax0 = 16.0
 vmax0 = 1.0
 xoffset = 7
	
 cutx  = 100
 dimx  = 400 - 2*cutx
 cutv  = 0
 dimv  = dimv0 - 2*cutv

 file1     = "phasex-1-sp0-3.250"
 file2     = "phasex-1-sp0-6.500"

 psfile    = "phase.ps"
 epsfile   = "phase.eps"

; ------ axes ----------------------------------------------------------------------------

 s1     = dimx
 s2     = dimv
 xmin   = xmax0/dimx0 * cutx - xoffset
 xmax   = xmin + xmax0/dimx0 * dimx 
 ymax   = vmax0 - vmax0/dimv0 * cutv  
 x      = xmin + (xmax-xmin)/s1 * findgen(s1) 
 y      = -ymax + 2.0*ymax/s2 * findgen(s2)

; ------ color table ---------------------------------------------------------------------

 max  = 255
 expo =	5
 r = findgen(max+1)
 g = findgen(max+1)
 for i=0, max do begin
    g(i) = max * ( r(i)/max ) ^ expo 
 endfor
 g(0) = 255
 g(1) = 0
 tvlct, g, g, g

; ======== device x ======================================================================

 set_plot,"x"
 window, 1, xpos=100, ypos=28, xsize=2*dimx+300, ysize=dimv+180, title="phasespace(x,v)"

; ------ image size ----------------------------------------------------------------------

 xoff = !d.x_px_cm*4 
 yoff = !d.y_px_cm*3
 s1x  = !d.x_px_cm*5
 s1y  = !d.y_px_cm*10
 s2x  = s1x    
 s2y  = s1y
 dist = s1x

; ---- load files ------------------------------------------------------------------------

 print, file1
 openu,1,file1
 a1 = assoc(1,bytarr(dimx0,dimv0,/nozero))
 phase1 = extrac(a1(0),cutx,cutv,dimx,dimv)

 print, file2
 openu,2,file2
 a2 = assoc(2,bytarr(dimx0,dimv0,/nozero))
 phase2 = extrac(a2(0),cutx,cutv,dimx,dimv)

; ---- plot ------------------------------------------------------------------------------

 tvscl, alog10(1+phase1), xsize=s1x, ysize=s1y, xoff, yoff

 tvscl, alog10(1+phase2), xsize=s2x, ysize=s2y, xoff+dist, yoff

 plot, x, y, xstyle=1, ystyle=1, xrange=[xmin,xmax], yrange=[-ymax,ymax],$
	xtitle="!3x/!7k!D0!N!3", ytitle="!3v!Dx!N/c!3",xticklen=+0.02,$
        yticklen=0.02, xticks=3, xtickv=[-2,0,2,4],$
	charsize=2, color=1, position=[xoff,yoff,xoff+s1x,yoff+s1y],$
	/noerase, /nodata, /device

 plot, x, y, xstyle=1, ystyle=1, xrange=[xmin,xmax], yrange=[-ymax,ymax],$
	xtitle="!3x/!7k!D0!N!3", ytitle="",xticklen=+0.02,$
        yticklen=0.02, xticks=3, xtickv=[-2,0,2,4],$
        yticks=4, ytickname=[' ',' ',' ',' ',' '],$
	charsize=2, color=1, position=[xoff+dist,yoff,xoff+dist+s2x,yoff+s2y],$
	/noerase, /nodata, /device

; ====== device ps =======================================================================

 set_plot,"ps"
 device, file=psfile, /color, bits=8
;device, /encapsulated, file=epsfile, /color, bits=8

; ------ image size ----------------------------------------------------------------------

 xoff = !d.x_px_cm*4 
 yoff = 0
 s1x  = !d.x_px_cm*5
 s1y  = !d.y_px_cm*10
 s2x  = s1x    
 s2y  = s1y
 dist = s1x

; ---- plot ------------------------------------------------------------------------------

 tvscl, alog10(1+phase1), xsize=s1x, ysize=s1y, xoff, yoff

 tvscl, alog10(1+phase2), xsize=s2x, ysize=s2y, xoff+dist, yoff

 plot, x, y, xstyle=1, ystyle=1, xrange=[xmin,xmax], yrange=[-ymax,ymax],$
	xtitle="!3x/!7k!D0!N!3", ytitle="!3v!Dx!N/c!3",xticklen=+0.02,$
        yticklen=0.02, xticks=3, xtickv=[-2,0,2,4],$
	charsize=2, color=1, position=[xoff,yoff,xoff+s1x,yoff+s1y],$
	/noerase, /nodata, /device

 plot, x, y, xstyle=1, ystyle=1, xrange=[xmin,xmax], yrange=[-ymax,ymax],$
	xtitle="!3x/!7k!D0!N!3", ytitle="",xticklen=+0.02,$
        yticklen=0.02, xticks=3, xtickv=[-2,0,2,4],$
        yticks=4, ytickname=[' ',' ',' ',' ',' '],$
	charsize=2, color=1, position=[xoff+dist,yoff,xoff+dist+s2x,yoff+s2y],$
	/noerase, /nodata, /device

 close,1
 close,2

 device, /close

 print, "picture written to file "+psfile 

; ----------------------------------------------------------------------------------------

end










