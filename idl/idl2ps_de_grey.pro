pro idl2ps_de

 yname   = "!3n!De!N!3/n!Dc!3"

 xoffset = 7.0000

 xmin    = 6.00 - xoffset
 xmax    = 10.00 - xoffset
 ymin    = 0
 ymax    = 6
 zmax    = 30

 level1  = 2.00e+02
 level2  = 3.00e+02
 level3  = 4.00e+02

 dim     = 400

 color   = 0

 unit    = "spacetime-de"

 idl2eps = "density-double-0-6.eps"

 openu,1,unit
 a=assoc(1,bytarr(dim,dim,/nozero))
 dens=a(0)
 openu,2,"scale-idl"
 b=assoc(2,bytarr(30,256,/nozero))
 sca=b(0)
 d=size(dens) 
 s=size(sca) 

 if ( color eq 0 ) then begin
	 for i=0, dim-1 do begin
		for j=0, dim-1 do begin
			if (dens(i,j) eq 0) then dens(i,j)=127
		endfor	
	 endfor		
	 max  = 255
	 max1 = 200
	 max2 = 128
	 expo1 = 1
	 expo2 = 1
	 r = findgen(max+1)
	 g = findgen(max+1)
	 for i=0, max2 do begin
	     g(i) = max * ( r(i)/max2 ) ^ expo1 
	 endfor
	 for i=max2+1, max1 do begin
	     g(i) = 255
	 endfor
	 for i=max1+1, max do begin
	     g(i) = max * ( 1 - (r(i)-max1)/(max-max1) ) ^ expo2 
	 endfor
	 g(0) = 255
	 tvlct, g, g, g
	 black = 255
	 white = 128
	 ccol_1  = black
	 ccol_2  = black
	 ccol_3  = black
 endif else begin
	 loadct,color
	 max  = 255
	 r = findgen(max+1)
	 g = findgen(max+1)
	 b = findgen(max+1)

	 tvlct, r, g, b, /get
 	 r(0) = 255
	 g(0) = 255
	 b(0) = 255
	 tvlct, r, g, b

	 black = 1
	 white = 0
	 ccol_1 = black
	 ccol_2 = black
	 ccol_3 = black
 endelse

 set_plot,"ps"
 device, /encapsulated, file=idl2eps, /color, bits=8
 xoff = !d.x_px_cm*2 
 yoff = 0
 s1x  = !d.x_px_cm*10
 s1y  = !d.y_px_cm*10
; s2x  = !d.x_px_cm    
 s2x  = !d.x_px_cm*0.5    
; s2y  = !d.y_px_cm*10
 s2y  = !d.y_px_cm*8
; dist = !d.x_px_cm*4
 dist = !d.x_px_cm

 tv, dens, xsize=s1x, ysize=s1y, xoff, yoff
 x = xmin + (xmax-xmin)/d(1)*findgen(d(1))
 y = ymin + (ymax-ymin)/d(2)*findgen(d(2))
 zcrit = 256/zmax
 contour, dens, /noerase, xstyle=4, ystyle=4, $
	  level = [level1*zcrit, level2*zcrit, level3*zcrit], $
          c_color = [ccol_1, ccol_2, ccol_3], $
          position=[xoff,yoff,xoff+s1x,yoff+s1y],/device
  plot, x, y, xstyle=1, ystyle=1, xrange=[xmin,xmax], yrange=[ymin,ymax],$
        xtitle="!3x/!7k!3!D0!3", ytitle="!3t/!7s!3",$
        charsize=1.75, color=black, /nodata, /noerase, /device,$
        position=[xoff,yoff,xoff+s1x,yoff+s1y]

; xoff = xoff + s1x + dist
; yoff = yoff
 xoff = xoff + !d.x_px_cm * 0.5
 yoff = yoff + !d.y_px_cm

; tv, sca, xsize=s2x, ysize=s2y, xoff, yoff
; x = 1.0/s(1)*findgen(s(1))
; y = zmax/s(2)*findgen(s(2))

; contour, sca, /noerase, xstyle=4, ystyle=1, yrange=[0,zmax], charsize = 2, $ 

; contour, sca, /noerase, xstyle=4, ystyle=4, yrange=[0,zmax], $ 
;    	  levels = [level1, level2, level3], $
;          c_color = [ccol_1, ccol_2, ccol_3], $
;	  color = black, $
;          position = [xoff,yoff,xoff+s2x,yoff+s2y],/device

; plot, x, y, xstyle = 4, ystyle = 1, yrange = [0,zmax],$
;          ytitle = yname, charsize = 2,$
;          color = black, /nodata, /noerase, /device,$
;          position = [xoff,yoff,xoff+s2x,yoff+s2y]

; axis, xaxis=0, color=black, xticks=1, xtickname=replicate(' ',2), xminor=1
; axis, xaxis=1, color=black, xticks=1, xtickname=replicate(' ',2), xminor=1
; axis, yaxis=0, color=black, yticks=1, ytickname=replicate(' ',2)
; axis, yaxis=1, color=black, yrange=[0,zmax], ystyle=1, ytitle=yname, ycharsize=1.5

 device, /close
 close,1
 close,2

end