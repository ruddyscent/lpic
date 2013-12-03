pro idl_inverted_grey_table
	
 max  = 255

 r = findgen(max+1)
 g = findgen(max+1)

 for i=1, max do begin
     g(i) = max * ( 1 - r(i)/max )  
 endfor
 
 g(0) = 255

 tvlct, g, g, g

end