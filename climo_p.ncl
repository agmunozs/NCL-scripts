;***************************************************************************************
;Centro de Modelado Científico (CMC). La Universidad del Zulia   -- cmc.org.ve
;Observatorio Andino de Eventos Extraordinarios (OA)             -- ole2.org
;Copyleft: Ángel G. Muñoz S. -- agmunoz@cmc.org.ve
;Descripción: script NCL del OA y CMC para generar archivo netcdf con climatologia del pais
;***************************************************************************************
;---------------------------------------------------------------------------------------
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"
;--------------------------------------------------------------------------------------- 

begin
  npts = 24 ; Num estaciones
  ntim = 324; Num meses  total 456. para la climo: 360
  minlat=5 ; lo obvio
  maxlat=11
  minlon=-86
  maxlon=-75
  reso=30  ; Resolucion deseada en km


;*****************************************************************************
;Seccion Lectura Archivo ASCII
;*****************************************************************************
  data = asciiread("data_PAN_P.dat",(/ntim*npts,6/),"float")

; Años y meses
  ano1d= data(::npts,0)  ;con el ::ALGO hacemos que lea cada ALGO pasos
  mes1d= data(::npts,1)

; Para las latidudes de las estaciones, basta adquirir las primeras 170 lineas:
  lat1d = data(0:npts-1,3)

; Idem para las longitudes
  lon1d  = data(0:npts-1,4)

;***************************************************************************
;Abrimos archivo CRU para first guess del analisis objetivo
;**************************************************************************
  f    = addfile ("precip_CRU_mes_71-00.nc", "r")
  pr   = f->prcp(35:,:,:)
  time = f->T(35:)
  lati = f->Y
  loni = f->X
  pr!0="time"
  pr!1="lat"
  pr!2="lon"
printVarSummary(time)

;*****************************************************************************
; Definimos la nueva malla para el Cressman
;*****************************************************************************
  lon = fspan(minlon, maxlon,floattoint(abs(minlon-maxlon)*111.111/reso)) ;Creamos la malla (longitud, res=30 km)
  nlon = dimsizes( lon)
  ;print(lon)
  lat = fspan(minlat, maxlat,floattoint(abs(maxlat-minlat)*111.111/reso))  ;Creamos la malla (latitud, res=30 km)
  nlat = dimsizes( lat)
 ; printVarSummary(lat)
 ; printVarSummary(lon)

;****************************************************
;Interpolacion a la malla del Cressman de la data CRU
;****************************************************
pre=linint2 (pr&lon, pr&lat, pr, False, lon, lat, 0)
pre@_FillValue=-99.9
pre@lon=lon
pre@lat=lat
pre!0="time"
pre!1="lat"
pre!2="lon"

printVarSummary(pre)
;print(pre)

;*****************************************************************************
; Definimos unos arreglos a usar
 temp1D=new((/ntim,npts/),typeof(data))
 temp2D=new((/ntim,npts,npts/),typeof(temp1D))

; Ciclo de lectura y almacenamiento en arreglos para la variable física 
  do i=0,ntim-1 
   temp1D(i,:)   = data(i*npts:(i+1)*npts-1,5)            
;   temp2D(i,:,:) = onedtond(temp1D(i,:),(/npts,npts/)) ; convertimos a 2D
  end do

;**************************************************
;Esto Define los valores faltantes (superfluo)
  temp1D@_FillValue    = -99.9


; Assign named dimensions 

 ; temp2D!0 = "lat"
 ; temp2D!1 = "lon"

; Assign coordinate variables

 ; temp2D&lat = lat1d
 ; temp2D&lon = lon1d

  ;printVarSummary(temp2D)
;============================================================== ;
  grid =new((/nlat,nlon,ntim/),typeof(temp1D))
  do nmo =0,ntim-1
    zVal = temp1D(nmo,:)
    rscan = (/0.5,0.3,.1/)
    opt = True
    opt@guess = pre(nmo,:,:)  ; CRU como fisrt guess, pasamos a mm
    ;opt@zonal = True ; esta opción (recomendada) permite que FG = medias zonales
    grid(:,:,nmo) = obj_anal_ic(lon1d,lat1d,zVal,lon,lat, rscan, opt) 
    ;printVarSummary(grid)
    ;print(grid)
  end do

; Ahora control de calidad mínima: valores negativos tras Cressman-->Dato faltante
  do it=0,ntim-1
   do i=0,nlat-1
     do j=0,nlon-1
       if ((.not.ismissing(grid(i,j,it))).and.(grid(i,j,it).lt.0)) then 
        grid(i,j,it)=-999
       end if
     end do
   end do
  end do
;**************************************************
; Atributos requeridos para evitar warnings fastidiosos
;**************************************************
  grid&lon@units    = "degrees_east"
  grid&lat@units     = "degrees_north"
  ;pre&lon@units    = "degrees_east"
  ;pre&lat@units     = "degrees_north"
;**************************************************
;Esto Define los valores faltantes (superfluo)
  ;grid@_FillValue    = 0.0
  grid@_FillValue    = -999


;******************************
;Ahora rellenamos con Poisson relajado
;********************************

  guess     = 1                ; use zonal means
  is_cyclic = False             ; cyclic [global]
  nscan     = 1500             ; usually much less than this
  eps       = 1.e-2            ; variable dependent
  relc      = 0.6              ; relaxation coefficient
  opt2       = 0                ; not used

 do nme=0,ntim-1
   poisson_grid_fill(grid(:,:,nme), is_cyclic, guess, nscan, eps, relc, opt2)
 end do
;**************************************************************************

;*********************************************************************************
;CLIMATOLOGIA y VARIABILIDAD

  Pclm = clmMonLLT (grid)       ; 
  Pstd = stdMonLLT (grid)       ;

;  Pclm&lat = lat
;  Pclm&lon = lon
;  Pstd&lat = lat
;  Pstd&lon = lon

  Pclm&lon@units    = "degrees_east"
  Pclm&lat@units    = "degrees_north"
  
  Pstd&lon@units    = "degrees_east"
  Pstd&lat@units    = "degrees_north"
  
  printVarSummary(Pclm)
  printVarSummary(Pstd)
  printVarSummary(Pclm&month)
;  months=Pclm&month
;  printVarSummary(months)


;******************
;Seccion de Graf
;******************
;*********************************************************************************

  wks = gsn_open_wks("ps","climo_p")             ; open a workstation
  wks2 = gsn_open_wks("ps","desviacion_p")
  
  gsn_define_colormap(wks,"WhBlGrYeRe");
  gsn_define_colormap(wks2,"WhBlGrYeRe")
;----------- Begin plot -----------------------------------------

  res                 = True       ; Indicate you want to set some
                                         ; resources.

  res@cnMonoLineColor = False      ; Turn off the drawing of
  res@gsnMaximize     = True
  
  res                          = True
  
  res@cnFillOn                 = True	             ; turn on color
  res@cnLinesOn                = False;True             ; turn off contour lines
  res@cnLevelSpacingF          = 20.              ; contour interval
  res@cnLevelSelectionMode     = "ManualLevels"
  res@cnMinLevelValF           = 10.
  res@cnMaxLevelValF           = 400.
  res@cnFillDrawOrder          = "Draw"         ; draw contours first
  ;res@lbAutoManage             = False
  ;res@lbLabelStride            = 3                 ; stride on label bar
  ;res@lbLabelFontHeightF       = 0.25
  res@gsnSpreadColors          = True             ; use full colormap
  res@gsnSpreadColorEnd        = -3                ; don't use land color
  res@gsnAddCyclic             = False             ; regional data 
  res@mpOutlineOn         = True
  res@mpMaxLatF           = maxlat
  res@mpMinLatF           = minlat
  res@mpMinLonF           = minlon
  res@mpMaxLonF           = maxlon
  res@gsnMajorLatSpacing  = 1
  res@gsnMajorLonSpacing  = 5
  res@gsnFrame            = False
  res@mpOutlineBoundarySets = "AllBoundaries"
  res@mpDataSetName = "Earth..4" 
  res@mpDataBaseVersion = "MediumRes"
  res@mpOceanFillColor            = "White"
;  res@mpFillAreaSpecifiers        = (/"Water"/)
;  res@mpSpecifiedFillColors       = (/"white"/)
  res@mpInlandWaterFillColor      = 238
  res@mpLandFillColor             = -1     ; transparent
;  res@mpMaskAreasSpecifiers       = (/"Water"/)
;  res@mpMaskOutlineSpecifiers     = (/"water"/)  ; mascara para agua

  txres               = True                     ; text mods desired
  txres@txFontHeightF = 0.01     


;************************************************
; create panel plots
;*************************************************
  plot   = new (4 , graphic)                    ; create graphical array

  res                      = True               ; plot options desired

  res@cnFillOn             = True               ; turn on color fill
  res@cnInfoLabelOn        = False              ; turn off contour info label
  res@cnLinesOn            = False              ; turn off contour lines
  res@cnLineLabelsOn       = False              ; turn off line labels
  ;res@cnLevelSelectionMode = "ExplicitLevels"   ; set explicit contour levels
  ;res@cnLevels             = (/ 0.5, 1.0, 2.0 \ ; set unequal contour levels
  ;                            , 3.0, 4.0, 5.0 \
  ;                            , 7.5,10.0 /)

  res@mpFillOn             = False              ; turn off gray continents
  ;res@mpCenterLonF         = 180                ; Centers the plot at 180 
 
  res@lbLabelBarOn         = False              ; No single label bar
  res@gsnDraw              = False
  res@gsnFrame             = False
  

  resP                     = True               ; panel options
  resP@txString            = "Climatologia ETESA - Analisis Cressman Precipitacion"   ; common title
  resP@gsnMaximize         = True               ; maximize image
  resP@gsnPanelLabelBar    = True               ; Add common label bar 
  resP@lbAutoManage             = False
  resP@lbLabelStride            = 3                 ; stride on label bar
  resP@lbLabelFontHeightF       = 0.015


month = (/"Enero", "Febrero", "Marzo", "Abril" \
          ,"Mayo", "Junio", "Julio", "Agosto"         \
          ,"Septiembre", "Octubre", "Noviembre"      \
          ,"Diciembre" /)

  i = -1                                        ; Climatologies
  do nmo=0,3                                 ; loop over the months
     i = i+1
     res@gsnCenterString   = month(nmo)+" 1974-2000";+time(0)/100 +"-"+ time(ntim-1)/100
     plot(i) = gsn_csm_contour_map_ce(wks,Pclm(:,:,nmo), res)  ; create plot
  end do

  gsn_panel(wks,plot,(/2,2/),resP)

 i = -1                                        ; Climatologies
  do nmo=4,7                                 ; loop over the months
     i = i+1
     res@gsnCenterString   = month(nmo)+" 1974-2000";+time(0)/100 +"-"+ time(ntim-1)/100
     plot(i) = gsn_csm_contour_map_ce(wks,Pclm(:,:,nmo), res)  ; create plot
  end do

  gsn_panel(wks,plot,(/2,2/),resP)

 i = -1                                        ; Climatologies
  do nmo=8,11                                 ; loop over the months
     i = i+1
     res@gsnCenterString   = month(nmo)+" 1974-2000";+time(0)/100 +"-"+ time(ntim-1)/100
     plot(i) = gsn_csm_contour_map_ce(wks,Pclm(:,:,nmo), res)  ; create plot
  end do

  gsn_panel(wks,plot,(/2,2/),resP)

;************************************************
; create second panel plot
;*************************************************

 resP@txString            = "Desviacion Estandard ETESA - Analisis Cressman Precipitacion"   ; common title

  i = -1                                        ; Std Deviations
  do nmo=0,3                                 ; loop over the months
     i = i+1
     res@gsnCenterString   = month(nmo)+" 1974-2000";+time(0)/100 +"-"+ time(ntim-1)/100
     plot(i) = gsn_csm_contour_map_ce(wks2,Pstd(:,:,nmo), res)  ; create plot
  end do

  gsn_panel(wks2,plot,(/2,2/),resP)


  i = -1                                        ; Std Deviations
  do nmo=4,7                                 ; loop over the months
     i = i+1
     res@gsnCenterString   = month(nmo)+" 1974-2000";+time(0)/100 +"-"+ time(ntim-1)/100
     plot(i) = gsn_csm_contour_map_ce(wks2,Pstd(:,:,nmo), res)  ; create plot
  end do

  gsn_panel(wks2,plot,(/2,2/),resP)

  i = -1                                        ; Std Deviations
  do nmo=8,11                                 ; loop over the months
     i = i+1
     res@gsnCenterString   = month(nmo)+" 1974-2000";+time(0)/100 +"-"+ time(ntim-1)/100
     plot(i) = gsn_csm_contour_map_ce(wks2,Pstd(:,:,nmo), res)  ; create plot
  end do

  gsn_panel(wks2,plot,(/2,2/),resP)


;--------------------------------------------------

;  do nme=0,ntim-1
;    xwks = gsn_open_wks("ps","ecuador_"+mes1d(nme)+"_"+ano1d(nme))
;    xwks@wkWidth = 600 
;    xwks@wkHeight= 800
;    gsn_define_colormap(xwks,"WhBlGrYeRe");
;    i = NhlNewColor(xwks,.96,.87,.7)                   ; land color
;    res@tiMainString    = "Precipitacion Acumulada (Analisis)  " +mes1d(nme)+"-"+ano1d(nme)
;    plot = gsn_csm_contour_map_ce(xwks,grid(:,:,nme),res)       ; Draw a contour plot.
;    ; Texto OA
;   gsn_text_ndc(xwks,"Climatologia con Analisis Objetivo: Cressman  --  First Guess: CRU",0.5,0.2,txres)
;   gsn_text_ndc(xwks,"Base de Datos: INAMHI -- Periodo: 1971-2000",0.5,0.175,txres)
;   gsn_text_ndc(xwks,"No. Estaciones: 169",0.5,0.15,txres)
;    gsn_text_ndc(xwks,"Centro de Modelado Cientifico (CMC) de La Universidad del Zulia -- Observatorio Andino (OA) -- http://cmc.org.ve -- INAMHI -- http://inamhi.gov.ec",0.5,0.1,txres) 
;    frame(xwks)                                     ; now advance frame
;  end do


 ;plot = gsn_contour(xwks,grid,res)  ; Draw a contour plot.


end
;------------------------------------------------------------------------------------- 
