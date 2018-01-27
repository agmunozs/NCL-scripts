load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"  
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"  
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl" 
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"

begin
;-------------------------------------------------------------
; MODIFICABLE POR USUARIO
;-------------------------------------------------------------
   ntim    = 970                     		     ;numero de tiempos (lineas en el archivo)
   fname   = "MaracayMayoJunio2011completo.txt"      ;Archivo de datos a leer (ASCII)

   estadf  = "estadisticoJunio.dat"        ;Archivo de descriptores estadisticos (SALIDA)
   cabef   = "cabeceraestaJunio.txt"       ;Archivo con cabecera de descriptores (SALIDA)
   seriesf = "icaroJunio"		     ;Archivo grafico para series de tiempo (SALIDA)
   boxf    = "boxplotJunio"		     ;Archivo grafico para boxplot (SALIDA)
   
   PLOT    = True                                    ;Si "False" solo escribe en ascii estadisticos     

   vName   = "Base Sucre vs ICARO. Mayo-Junio  2011"      ;Titulo grafica serie de tiempo
   omm     = "SUCRE"                                 ;Nombre de la estacion OMM

   extgraf = "png"			             ;Extension del archivo grafico (e.g. ps, pdf, ...)
   niveles = 5                                       ;Niveles a dividir el intervalo en las figuras de sesgo
   ;Siguen los colores para los marcadores de las figuras de sesgo
   colors  = (/"white","black","firebrick","red","blue","blue4"/)

;*****************************************************************************
;FIN INTERVENCION USUARIO
;*****************************************************************************
   if (PLOT) then
       wksType = extgraf
       wksName = seriesf   
   end if

  obs     = new((/3,ntim/),float)  ;estas dos tienen 3 indices, uno para T,HR,P
  ica     = new((/3,ntim/),float)
  delticas= new((/3/),float)
;*********************************
;LECTURA
;*********************************
;Si hiciera falta, modificar los numeros al final de cada str_get_cols. Son columnas de inicio y fin para cada variable
  data    = asciiread(fname,-1,"string")
  year    = stringtointeger(str_get_cols(data,6,9))
  month   = stringtoint(str_get_cols(data,0,1))
  day     = stringtoint(str_get_cols(data,3,4))
  hr      = stringtoint(str_get_cols(data,11,12))
  mi      = stringtoint(str_get_cols(data,14,15))
  obs(0,:)= stringtofloat(str_get_cols(data,17,21))
  obs(1,:)= stringtofloat(str_get_cols(data,29,32))
  obs(2,:)= stringtofloat(str_get_cols(data,37,42))
  ica(0,:)= stringtofloat(str_get_cols(data,46,50))
  ica(1,:)= stringtofloat(str_get_cols(data,64,69))
  ica(2,:)= stringtofloat(str_get_cols(data,78,86))


;**********************************
;VARIABLES DE TIEMPO
;**********************************
   ddd     = day_of_year(year, month, day)
   dfrac   = hr 
   yyyyddd = year*1000000 + month*10000 +day*100+hr                  
   tempo   = yyyymmddhh_to_yyyyfrac(yyyyddd,0)
   ;print(tempo+"   "+yyyyddd)
   hhmm    = tempo;day+hr/24.+mi/86400.

;**********************************
;ESTADISTICOS
;**********************************
  opt   = False
  do p=0,2
   delticas(p)=dim_avg(ica(p,:))-dim_avg(obs(p,:))
  end do
  o1    = stat_dispersion(obs(0,:), opt )
  i1    = stat_dispersion(ica(0,:)-delticas(0), opt )
  o2    = stat_dispersion(obs(1,:), opt )
  i2    = stat_dispersion(ica(1,:)-delticas(1), opt )
  o3    = stat_dispersion(obs(2,:), opt )
  i3    = stat_dispersion(ica(2,:)-delticas(2), opt )
  delta1=i1-o1
  delta2=i2-o2
  delta3=i3-o3

  nStat = num(.not.ismissing(o1))
  statLabel = new ( nStat, "string")

  statLabel = (/ "   Media", "   sigma", "     Min" \
               , "Dec.infe", "Oct.infe", "Sex.infe" \
               , "Cuartinf", "Tri.infe" \
               , " Mediana" \
               , "Tri.supe", "Cuartsup", "Sex.supe","Oct.supe", "Dec.supe" \
               , "     Max", "   Rango","  Disper","     RMS" \
               , "    nTot", "    nUso","   nFalt","   %Falt" /)

  encabezado="Descrip.       OBS:T       ICA:T      DELTA      OBS:HR      ICA:HR      DELTA      OBS:P       ICA:P      DELTA"
  salida=(/statLabel+"   " \
       +sprintf("%9.2f", o1(0:nStat-1))+"   " \
       +sprintf("%9.2f", i1(0:nStat-1))+"   " \
       +sprintf("%9.2f", delta1(0:nStat-1))+"   " \
       +sprintf("%9.2f", o2(0:nStat-1))+"   " \
       +sprintf("%9.2f", i2(0:nStat-1))+"   " \
       +sprintf("%9.2f", delta2(0:nStat-1))+"   " \
       +sprintf("%9.2f", o3(0:nStat-1))+"   " \
       +sprintf("%9.2f", i3(0:nStat-1))+"   " \
       +sprintf("%9.2f", delta3(0:nStat-1))/)

  asciiwrite(estadf,salida)
  asciiwrite(cabef,encabezado)


;************************************************
; GRAFICOS 
;*******************************************G****
   if (PLOT) then
;*****************
;SERIES TEMPORALES
;*****************
        plot  = new ( 6, graphic)
        data2  = new ( (/2,ntim/), float, -999)

        wks   = gsn_open_wks (wksType,wksName)    

        res                   = True                    
        res@gsnDraw           = False
        res@gsnFrame          = False
        res@trXMinF           = min(hhmm)
        res@trXMaxF           = max(hhmm)
        ;res@trYMinF           = 0
        ;res@tiMainString      = ""
       
        res@xyLineThicknesses = (/1.0, 1.0/)              
        res@xyLineColors      = (/"blue","firebrick"/)          
        res@xyMonoDashPattern = True                      

        res@xyMarkLineMode    = "MarkLines"                
        res@xyMonoMarkerColor = False
        res@xyMarkers         = (/16,6/)                      
        res@xyMarkerColors    = (/"blue","firebrick"/)           
        res@xyMarkerSizes     = (/0.01,0.01/)                   
        ;res@tmXBPrecision     = 2
        res@gsnStringFontHeightF=0.03
	uni=(/"~S~o~N~C","%","mb"/)
        tit=(/"Temperatura","Humedad Relativa","Presion Superficial","~F33~D~F21~T","~F33~D~F21~HR","~F33~D~F21~P"/)

;print(num(ismissing(obs(0,:))))

;VARIABLES
        res@tiXAxisString         = "Fecha (HLV)"
        res@tmXBMode = "Explicit"	
        res@tmXBValues = (/2011.410,2011.420,2011.430,2011.440/)
        res@tmXBLabels = (/"28 May","01 Jun","05 Jun","09 Jun"/)
        res@tmXBMinorValues  = (/2011.405,2011.4075,2011.4125,2011.415,2011.4175,2011.425,2011.415,2011.4175,2011.4225,2011.425,2011.4275,2011.4325,2011.435,2011.4375,2011.445/)


        do i=0,2 
           deltica=dim_avg(ica(i,:))-dim_avg(obs(i,:))
           data2(0,:) = obs(i,:)
           data2(1,:) = ica(i,:)-deltica

           res@gsnCenterString       = tit(i)
           res@tiYAxisString         = uni(i)
           if(i.eq.2) then
             res@trYReverse=True
 	   end if
  	   if(i.eq.1) then
 	     res@trYMaxF           = 100
	   end if
           plot(i)  = gsn_csm_xy (wks,hhmm,data2,res) 
	   delete(res@trYMaxF)
        end do
        delete(res@trYReverse)
	delete(res@xyMarkers)
 	delete(res@xyMarkerColors)
;SESGOS
        res@xyMarkLineMode         = "Markers"
        res@xyMonoMarkerColor      = False
        res@xyMarkerThicknessF     = 2.5
	res@xyMarkerColors         = colors(2:)  
	res@gsnYRefLine            = 0.0            
        res@pmLegendDisplayMode    = "Always"         ; Turn on the legend
	res@lgAutoManage	   = False
        res@lgOrientation          = "horizontal"     ; Default is vertical
        res@pmLegendWidthF         = 0.70             ; Make it wider
        res@pmLegendOrthogonalPosF = -0.1             ; Move it up slightly
        res@lgPerimOn              = False            ; Turn off the perimeter box
	res@xyMarkers              = ispan(2,16,1)

        do j=3,5
           sesgo=(ica(j-3,:)-delticas(j-3))-obs(j-3,:)
 	   
           mnmxint = nice_mnmxintvl(min(sesgo),max(sesgo),niveles,True)
           nlevels = toint(((mnmxint(1)-mnmxint(0))/mnmxint(2))) + 1
           levels  = fspan(mnmxint(0),mnmxint(1),nlevels)
           npts=ntim;-num(ismissing(sesgo)) 
 	   data2d = new((/nlevels-1,npts/),typeof(sesgo))

       	   ;---Group the values and put in 2D array.
  	   labels = new(nlevels-1,string)
  	   do i=0,nlevels-2
    	     ii = ind(sesgo.ge.levels(i).and.sesgo.lt.levels(i+1))
             data2d(i,ii) = sesgo(ii)
	     labels(i) = levels(i) + ":" + levels(i+1)
    	     delete(ii)
  	   end do
	   delete(levels)

  	   res@xyExplicitLabels      = labels 
    	   res@lgLabelFontHeightF    = .03  
	   res@lgLabelJust	     = "CenterCenter"      
           res@gsnCenterString       = tit(j)
           res@tiYAxisString         = uni(j-3)
           plot(j)  = gsn_csm_xy (wks,hhmm,data2d,res)
	   delete(data2d)
           delete(labels)
           delete(res@xyExplicitLabels)
 	end do


  ;print("min/max data = " + min(data1d) + "/" + max(data1d))
  ;print("# of groups  = " + (nlevels-1))
  ;print("levels start at = " + levels(0) + ", end at " + levels(nlevels-1) + \
  ;      ", with a spacing of " + mnmxint(2))


        resP                     = True                
        resP@txString            = vName
        resP@gsnMaximize         = True
        resP@gsnPaperOrientation = "landscape"

        gsn_panel(wks,plot,(/2,3/),resP)              

;**********
;BOXPLOTS
;**********

  boxT     = new((/2,5/),"float",-999)
  boxHR    = new((/2,5/),"float",-999)
  boxP     = new((/2,5/),"float",-999)

  boxT(0,0)=o1(2)          	; min
  boxT(0,1)=o1(6)		; 1er-cuatil	
  boxT(0,2)=o1(8)		; mediana
  boxT(0,3)=o1(10)		; 3er cuartil
  boxT(0,4)=o1(14)		; max
  boxT(1,0)=i1(2)
  boxT(1,1)=i1(6)
  boxT(1,2)=i1(8)
  boxT(1,3)=i1(10)
  boxT(1,4)=i1(14)

  boxHR(0,0)=o2(2)
  boxHR(0,1)=o2(6)
  boxHR(0,2)=o2(8)
  boxHR(0,3)=o2(10)
  boxHR(0,4)=o2(14)
  boxHR(1,0)=i2(2)
  boxHR(1,1)=i2(6)
  boxHR(1,2)=i2(8)
  boxHR(1,3)=i2(10)
  boxHR(1,4)=i2(14)

  boxP(0,0)=o3(2)
  boxP(0,1)=o3(6)
  boxP(0,2)=o3(8)
  boxP(0,3)=o3(10)
  boxP(0,4)=o3(14)
  boxP(1,0)=i3(2)
  boxP(1,1)=i3(6)
  boxP(1,2)=i3(8)
  boxP(1,3)=i3(10)
  boxP(1,4)=i3(14)

  x=(/1,2/)
  wks2 = gsn_open_wks(extgraf,boxf)
  res2            = True                         ; plot mods desired
  res2@tmXBLabels = (/omm+":T","ICA:T"/) ; labels for each box
  res2@tiYAxisString         = uni(0)
  ;res2@tiMainString = "BOXPLOT: "+vName
  llres                   = True			
  llres@gsLineThicknessF  = 2.5                 ; line thickness 
  opti          = True			
  opti@boxWidth = .25				; Width of box (x units)
  opti@boxColors = (/"blue","red"/)  	; Color of box(es)
  plot1 = boxplot(wks2,x,boxT,opti,res2,llres)	; All 3 options used...
  res2@tmXBLabels = (/omm+":HR","ICA:HR"/) ; labels for each box
  ;res2@tiMainString = "BOXPLOT: "+vName
  res2@tiYAxisString         = uni(1)
  plot2 = boxplot(wks2,x,boxHR,opti,res2,llres)   ;
  res2@tmXBLabels = (/omm+":P","ICA:P"/) ; labels for each box
  ;res2@tiMainString = "BOXPLOT: "+vName
  res2@tiYAxisString         = uni(2)
  plot3 = boxplot(wks2,x,boxP,opti,res2,llres)   ;
  resP                     = True
  resP@txString            = "BOXPLOT: "+vName
  resP@gsnMaximize         = True
  resP@gsnPaperOrientation = "landscape"

  mres               = True                     ; marker mods desired
  mres@gsMarkerIndex = 3                        ; polymarker style
  mres@gsMarkerSizeF = 20.                      ; polymarker size
  mres@gsMarkerColor = "blue"                    ; polymarker color
  xi  = o1(0)                               ; polymarker locations
  xi2 = o2(0)
  xi3 = o3(0)

  dum1 = gsn_add_polymarker(wks2,plot1,1.,xi,mres) 
  dum2 = gsn_add_polymarker(wks2,plot2,1.,xi2,mres) 
  dum3 = gsn_add_polymarker(wks2,plot3,1.,xi3,mres) 
  
  mres@gsMarkerColor = "red"              ; change color
  mres@gsMarkerIndex = 3                        ; change style
  mres@gsMarkerSizeF = 20.                      ; change size
  xi  = i1(0)                               ; polymarker locations
  xi2 = i2(0)
  xi3 = i3(0)

  dum4 = gsn_add_polymarker(wks2,plot1,2.,xi,mres) 
  dum5 = gsn_add_polymarker(wks2,plot2,2.,xi2,mres)
  dum6 = gsn_add_polymarker(wks2,plot3,2.,xi3,mres) 
  gsn_panel(wks2,(/plot1,plot2,plot3/),(/1,3/),resP)

;**********
;HISTOGRAMAS
;**********
  wks3  = gsn_open_wks (extgraf,"histo")
  res3                         = True
  res3@gsnHistogramHorizontal  = True            ; draw horizontal histogram
  res3@gsnHistogramNumberOfBins = 10
  res3@gsnHistogramComputePercentagesNoMissing=True
  res3@gsnHistogramPercentSign= True
  res3@gsFillColor = colors(2:)
  res3@gsnDraw           = False
  res3@gsnFrame          = False

  res3@tiYAxisString         = "~F33~D~F21~T ("+uni(0)+")"
  res3@tiXAxisString	     = "Frecuencia Absoluta"
  plot4 = gsn_histogram(wks3,ica(0,:)-delticas(0)-obs(0,:),res3) 

  res3@tiYAxisString         = "~F33~D~F21~HR ("+uni(1)+")"
  plot5 = gsn_histogram(wks3,ica(1,:)-delticas(1)-obs(1,:),res3)

  res3@tiYAxisString         = "~F33~D~F21~P ("+uni(2)+")"
  plot6 = gsn_histogram(wks3,ica(2,:)-delticas(2)-obs(2,:),res3)

  resP@txString            = "HISTOGRAMAS: "+vName
  gsn_panel(wks3,(/plot4,plot5,plot6/),(/1,3/),resP)


end if

k   = ind( .not.(ismissing(obs(0,:)) .or. ismissing(ica(0,:))) )
  spc = spcorr( obs(0,k), ica(0,k))  
  signif=(1-rtest(spc, dimsizes(k), 0))
print("Correlacion Sp Temperatura:"+spc+" Significancia (alfa)="+signif)
k   = ind( .not.(ismissing(obs(1,:)) .or. ismissing(ica(1,:))) )
  spc = spcorr( obs(1,k), ica(1,k))  
  signif=(1-rtest(spc, dimsizes(k), 0))
print("Correlacion Sp Humedad:"+spc+" Significancia (alfa)="+signif)
kk   = ind( .not.(ismissing(obs(2,:)) .or. ismissing(ica(2,:))) )
  spc = spcorr( obs(2,kk), ica(2,kk))  
  signif=(1-rtest(spc, dimsizes(kk), 0))
print("Correlacion Sp Presion Superficial:"+spc+" Significancia (alfa)="+signif)
end
