xmax=7
xmin=12
ymin=0.000001
ymax=1
xtitle="log!d10!n(M/M"
ytitle='!4U!3/(Mpc!u-3!nlog!d10!nM!u-1!n)'
;cgSymbol('Sun')

;cggreek('phi')
;sunsymbol(font=0)

set_plot, 'PS'
!p.charsize=2.9
!p.charthick = 9
!p.thick = 9
!x.thick = 9
!y.thick = 9
!z.thick = 9
xsize=40		;14 per graph
aspect_ratio=1.0
ysize = xsize/aspect_ratio
device, filename = 'smf_lowz.eps', xsize=xsize, ysize=ysize, /color, xoffset=1, yoffset=5, /encapsulated 
multiplot,/init
multiplot,[2,2]
;multiplot,[2,2],mYtitle=ytitle,mxtitle=xtitle+cgSymbol('Sun')+')',mxTitSize=2.0,myTitSize=2.0
  
  
 readcol,'./Henriques2013a_smf_z0.00.txt',x,y
 plot, x,(10^(y))*(0.703^3),xrange=[xmin,xmax], yrange=[ymin,ymax], linestyle =0, color=0,/ylog,ytitle=ytitle,$
 	xmargin=[0,14],ymargin=[0,14],charsize=2.9
 readcol,'./obs_z0.txt',mass,phi,error
 symbols, 2, 1.0
 oploterror, mass, phi*(0.703^3), error*(0.703^3), color = 4, errcolor = 4, psym = 8,HATLENGTH = 100.0
 xyouts, 11.0, 0.1, 'z=0',charsize=2.9
  legend,['Baldry 2008','Li & White 2009'], /bottom, /right, charsize=2.5

 multiplot
 
 
  readcol,'./Henriques2013a_smf_z1.00.txt',x,y
 plot, x,(10^(y))*(0.703^3),xrange=[xmin,xmax], yrange=[ymin,ymax], linestyle =0, color=0,/ylog,$
 	xmargin=[0,14],ymargin=[0,14],charsize=2.9
 readcol,'./obs_z1.txt',mass,phi,error

 oploterror, mass, phi*(0.703^3), error*(0.703^3), color = 4, errcolor = 4, psym = 8,HATLENGTH = 100.0
  xyouts, 11.0, 0.1, 'z=1',charsize=2.9
   legend,['Fontana 2006','Perez-Gonzalez 2008','Ilbert 2010','Pozzetti 2010'], /bottom, /right, charsize=2.5

 multiplot


 readcol,'./Henriques2013a_smf_z2.00.txt',x,y
 plot, x,(10^(y))*(0.703^3),xrange=[xmin,xmax], yrange=[ymin,ymax], linestyle =0, color=0,/ylog,xtitle=xtitle+cgSymbol('Sun')+')',ytitle=ytitle,$
 	xmargin=[0,14],ymargin=[0,14],charsize=2.9
 readcol,'./obs_z2.txt',mass,phi,error
 symbols, 2, 1.0
 oploterror, mass, phi*(0.703^3), error*(0.703^3), color = 4, errcolor = 4, psym = 8,HATLENGTH = 100.0
  xyouts, 11.0, 0.1, 'z=2',charsize=2.9
   legend,['Marchesini 2009','Ilbert 2010','Sanchez 2011'], /bottom, /right, charsize=2.5

 multiplot


 readcol,'./Henriques2013a_smf_z3.00.txt',x,y
 plot, x,(10^(y))*(0.703^3),xrange=[xmin,xmax], yrange=[ymin,ymax], linestyle =0, color=0,/ylog,xtitle=xtitle+cgSymbol('Sun')+')',$
 	xmargin=[0,14],ymargin=[0,14],charsize=2.9
 readcol,'./obs_z3.txt',mass,phi,error
 symbols, 2, 1.0
 oploterror, mass, phi*(0.703^3), error*(0.703^3), color = 4, errcolor = 4, psym = 8,HATLENGTH = 100.0
 xyouts, 11.0, 0.1, 'z=3',charsize=2.9
  legend,['Marchesini 2009','Marchesini 2010','Sanchez 2011'], /bottom, /right, charsize=2.5

 multiplot

 
 
 
 

device, /close_file
set_plot,'x'




