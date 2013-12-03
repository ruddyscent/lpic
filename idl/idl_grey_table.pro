pro idl_grey_table
	
 max  = 255
 max1 = 213
 max2 = 128

 r = findgen(max+1)
 g = findgen(max+1)

 for i=0, max2 do begin
     g(i) = max * ( r(i)/max2 )
 endfor

 for i=max2+1, max1 do begin
     g(i) = 255
 endfor

 for i=max1+1, max do begin
     g(i) = max * ( 1 - (r(i)-max1)/(max-max1) )  
 endfor

 g(0) = 255

 tvlct, g, g, g

end