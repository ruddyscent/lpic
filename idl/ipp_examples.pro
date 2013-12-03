pro ipp_init
  COMMON ex_com,ex_path,image1,example
    erase
    idlpath=getenv('IDL_DIR')
    ex_path=idlpath+'/examples/ipp/'
    read_gif,ex_path+'/images/rzg.gif',image1,r,g,b
    tvlct,r,b,g
    tv,image1,100,50
end
pro make_example_event,event
  COMMON ex_com
    widget_control,event.id,get_uvalue=evval
       case evval of
          'DONE' : begin
		      widget_control,event.top,/destroy
		      ipp_init
		   end
          'RUN'  : case example of
                      'for_idl' : begin
				     cmd='xterm -e '+ex_path+'for_idl/for_idl'
				     spawn,cmd
				  end
                      'c_idl'   : begin
				     cmd='xterm -e '+ex_path+'c_idl/c_idl'
				     spawn,cmd
				  end
                      else      : dumy=execute(example)
		   endcase
          'CODE' : begin
		     dummy=''
                     pfile=pickfile(/read,path=ex_path+example)
		     if pfile ne dummy then xdisplayfile,pfile
		   end
       endcase 
end
pro make_example
  COMMON ex_com
    ex_plotbase=widget_base(/column)
    XPdMenu,[ '"     Done     "    DONE' , $
              '"     RUN      "	   RUN'    $
	    ] ,ex_plotbase
    case  example of
       'for_idl' : button=widget_button(ex_plotbase,value='display code',uvalue='CODE')
       'c_idl'   : button=widget_button(ex_plotbase,value='display code',uvalue='CODE')
       else      : dummy=widget_label(ex_plotbase,value='   ')
    endcase
    widget_control,ex_plotbase,/realize
    xmanager,'make_example',ex_plotbase
end
pro ipp_examples_event,event
  COMMON ex_com
    widget_control,event.id,get_uvalue=evval
       case evval of
       'DONE'     : widget_control,event.top,/destroy
       'RESET'    : ipp_init
       'C_IDL'    : begin
		      erase
		      read_gif,ex_path+'/images/c_idl.gif',image2,r,g,b
 		      tvlct,r,b,g
  		      tv,image2
	              example='c_idl'
		      make_example
	            end
       'FOR_IDL'  : begin
		      erase
		      read_gif,ex_path+'/images/for_idl.gif',image2,r,g,b
 		      tvlct,r,b,g
  		      tv,image2
	              example='for_idl'
		      make_example
	            end
       'EX_PLOT'  : begin
		      erase
		      read_gif,ex_path+'/images/text_plot.gif',image2,r,g,b
 		      tvlct,r,b,g
  		      tv,image2
	              example='ex_plot'
		      make_example
	            end
       'IDL_FOR'  : begin
		      erase
		      read_gif,ex_path+'/images/idl_for.gif',image2,r,g,b
 		      tvlct,r,b,g
  		      tv,image2
	              example='idl_for'
		      make_example
	            end
       'HIS_PLOT' : begin
		      erase
		      read_gif,ex_path+'/images/his_plot.gif',image2,r,g,b
 		      tvlct,r,b,g
  		      tv,image2
	              example='his_plot'
		      make_example
	            end
       'CODE'     : XDisplayFile,ex_path+'ipp_examples.pro'
       endcase 

end
pro ipp_examples
  COMMON ex_com
    base=widget_base(column=2)
    base1=widget_base(base,/column)
    base2=widget_base(base,row=2)
    lab2=widget_label(base2,xsize=600,ysize=80,/align_center,$
	              value="IPP-Examples Overview")
    XPdMenu,[ '" Done "		DONE' , $
              '" Reset "	RESET'  $
	    ],base1
    lab1=widget_label(base1,xsize=80,ysize=40,value="  ")
    butt1=widget_button(base1,value='  Histogram Plot  ',uvalue='HIS_PLOT')
    label=widget_label(base1,ysize=3,value='   ')
    butt2=widget_button(base1,value='  Plot with Text  ',uvalue='EX_PLOT')
    label=widget_label(base1,ysize=3,value='   ')
    butt3=widget_button(base1,value='  IDL calls Fortran  ',uvalue='IDL_FOR')
    label=widget_label(base1,ysize=3,value='   ')
    butt4=widget_button(base1,value='  Fortran calls IDL  ',uvalue='FOR_IDL')
    label=widget_label(base1,ysize=3,value='   ')
    butt5=widget_button(base1,value='  C calls IDL  ' ,uvalue='C_IDL')
    label=widget_label(base1,ysize=3,value='   ')
    butt6=widget_button(base1,value=' display code ',uvalue='CODE')
    draw=widget_draw(base2,xsize=800,ysize=800,colors=256,retain=2)
    widget_control,base,/realize
    ipp_init
    xmanager,'ipp_examples',base
end
