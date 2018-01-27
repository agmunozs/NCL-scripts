; ===============================================
; CMC - Perfilador Vertical en función del tiempo para Temperatura
; Escrito para visualizar perfiles de Ícaro tras Expedición Catatumbo 2015             
; ===============================================
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"

begin
;================================================
; get data
;================================================
  ntim    = 129                     		 ;numero de tiempos (lineas en el archivo)
  fnames  = (/"Omega-Cata15_tarde_Camp1.txt","Omega-Cata15_tardeSUP_Camp1.txt"/)
  fname   = fnames(0)      ;Archivo de datos a leer (ASCII)

   estadf  = "estadisticos.dat"        ;Archivo de descriptores estadisticos (SALIDA)
   cabef   = "cabecera.txt"            ;Archivo con cabecera de descriptores (SALIDA)
   seriesf = "Camp1_tarde"		       ;Archivo grafico para series de tiempo (SALIDA)
   boxf    = "boxplot"		       ;Archivo grafico para boxplot (SALIDA)
   
   PLOT    = True                                    ;Si "False" solo escribe en ascii estadisticos     

   vName   = "Campamento 1. 13-04-2015"      ;Titulo grafica serie de tiempo
   omm     = "Campamento1"                                 ;Nombre de la estacion OMM

   extgraf = "png"
   nlev    = 2
;*****************************************************************************
;FIN INTERVENCION USUARIO
;*****************************************************************************
  err = NhlGetErrorObjectId()
  setvalues err
  "errLevel" : "Fatal"                                                  ;only report Fatal errors
  end setvalues
  
  if (PLOT) then
       wksType = extgraf
       wksName = seriesf   
  end if

  nlevf     = dimsizes(fnames)
  obs       = new((/nlevf,3,ntim/),float)  ;estas dos tienen 3 indices, uno para T,HR,P
  slice     = new((/nlev,ntim/),float)
  slice!0   = "lev"         ; assign named dimension so can reorder and plot
  slice!1   = "time"
  slice2     = slice
  tmpr      = new((/ntim,2/),float)
  T = tmpr
  RH = tmpr
;*********************************
;LECTURA
;*********************************
;Si hiciera falta, modificar los números al final de cada str_get_cols. Son columnas de inicio y fin para cada variable
  data    = asciiread(fname,-1,"string")
  year    = stringtointeger(str_get_cols(data,0,3))
  month   = stringtoint(str_get_cols(data,5,6))
  day     = stringtoint(str_get_cols(data,8,9))
  hr      = stringtoint(str_get_cols(data,11,12))
  mi      = stringtoint(str_get_cols(data,14,15))

delete(data)
do fn=0,nlevf-1
  data       = asciiread(fnames(fn),-1,"string")
  obs(fn,0,:)= stringtofloat(str_get_cols(data,37,40))
  obs(fn,1,:)= stringtofloat(str_get_cols(data,41,44))
  obs(fn,2,:)= stringtofloat(str_get_cols(data,46,51))*68.947 ;Omega mide en psi… transformando a mb
  ;if (fn.eq.1)
  ;  obs(fn,2,:)=1009.
  ;  obs(fn,0,:)=26.
  ;  obs(fn,1,:)=94.
  ;end if
  delete(data)
end do

do it = 0,ntim-1
   tmpr(it,:) = (/obs(0,2,it),obs(1,2,it)/)
   T   (it,:) = (/obs(0,0,it),obs(1,0,it)/)
   RH  (it,:) = (/obs(0,1,it),obs(1,1,it)/)
end do
;print(obs);
;print(tmpr)
print(max(tmpr))
print(min(tmpr))
;**********************************
;VARIABLES DE TIEMPO
;**********************************
   ddd     = day_of_year(year, month, day)
   dfrac   = hr 
   yyyyddd = year*1000000 + month*10000 +day*100+hr                  
   tempo   = yyyymmddhh_to_yyyyfrac(yyyyddd,0)
   ;print(tempo+"   "+yyyyddd)
   hhmm    = day+hr/24.+mi/86400.

;************************************************
; INTERPOLACIÓN VERTICAL
;************************************************
  lev        = fspan(1013.,980.,nlev);fspan(max(tmpr),min(tmpr),nlev)
  lev@units  = "mb"        ; units attribute required by gsn_csm_pres_hgt

do it = 0,ntim-1
  slice (:,it)= int2p(tmpr(it,:),T(it,:),lev,1)
  slice2(:,it)= int2p(tmpr(it,:),RH(it,:),lev,1)
end do

  slice&lev = lev           ; required by gsn_csm_pres_hgt
  slice2&lev = lev           ; required by gsn_csm_pres_hgt
  ;slice&time= tempo
 print(slice) 

;================================================
; smooth data
;================================================
  wgt   = (/ 1., 3., 4., 3., 1./)                  ; wgts for temporal smooth
  wgt   = wgt/sum(wgt)                             ; normalize

  slice   = wgt_runave(slice, wgt, 0) 
  slice2  = wgt_runave(slice2, wgt, 0) 
;================================================
; plot
;================================================
  wks   = gsn_open_wks (wksType, wksName)            ; open ps file
  gsn_define_colormap(wks,"BlRe")                ; choose colormap
  resL                      = True                  ; plot mods desired
  
  resL@cnFillOn             = True                  ; turn on color
  resL@cnLinesOn            = False                 ; no contour lines

  ;resL@cnLevelSelectionMode = "ManualLevels"       ; set manual contour levels
  ;resL@cnMinLevelValF       = 23.                  ; set min contour level
  ;resL@cnMaxLevelValF       = 32.                  ; set max contour level
  ;resL@cnLevelSpacingF      = 0.5                 ; set contour spacing
  
  resL@tiMainString = "Congo Mirador - 13/14.04.2015"
  resL@tiMainFontHeightF = 0.02
  resL@vpXF = 0.13 
  resL@vpWidthF = 0.70 
  resL@vpHeightF = 0.45

; in this example, we only plot part of the vertical extent. This reduces
; the number of tickmarks on the height axis. The following will increase
; the number:

  resL@tmYRTickSpacingF     = 2                   ;
       ;resL@tiXAxisString         = "Fecha (HLV)"

        resL@tmXBMode = "Explicit"	
        resL@tmXBValues = (/3,15,27,39,51,63,75,87,99,111,123/)
        resL@tmXBLabels = (/"19h","20h","21h","22h","23h","00h","01h","02h","03h","04h","05h"/)
        ;resL@tmXBValues = (/6,12,18,24,30/)
        ;resL@tmXBLabels = (/"11h","12h","13h","14h","16h"/)
        
      ;resL@tmYLMode = "Explicit"	
      ;resL@tmYLValues = (/1000.,990.,980.,970.,960.,950./)
      ;resL@tmYLLabels = (/“0”,”200","280","366”,”452”,”540”,”700”/)
      ;resL@tiYAxisString         = "Altura (m)"

        resL@tmYLMode = "Explicit"	
        resL@tmYLValues = (/1000.,990.,980.,970.,960.,950.,940.,930.,920./)
	resL@tmYLLabels = (/"1000","990","980","970","960","950","940","930","920"/)
        resL@tiYAxisString         = "Presio~H-13V2F35~B~FV-2H3~n (mb)"
   	resL@tiXAxisFontHeightF    = 0.015
	resL@lbTitleString = "Temperatura (Celsius, colores) y Humedad Relativa (%, contornos)"
	resL@lbTitleFontHeightF = 0.012
	resL@lbTitlePosition = "Bottom"
 	
  resL@gsnDraw              = False                 ; don't draw yet
  resL@gsnFrame             = False                 ; don't advance frame yet
  resL@cnInfoLabelOn = False

  plot = gsn_csm_pres_hgt (wks,slice(:,0:ntim-1),resL) ;


  ;resL@tmYUseLeft = True
  resL@cnFillOn             = False                  ; turn on color
  resL@cnLinesOn            = True
  ;    resL@tmYLMode = "Explicit"	
  ;    resL@tmYLValues = (/1000.,990.,980.,970.,960.,950./)
  ;    resL@tmYLLabels = (/“0”,”200","200","300","400","500","1200"/)
  ;    resL@tiYAxisString         = "Altura (m)"
  ; to remove the "height" label on the right, we have to go through several
; steps. The label is not an axis string but extra text placed there by
; the plot template
  

  getvalues plot@contour                                                     
  "pmAnnoManagers" : am_ids                                                
  end getvalues                                                              
  
  index = ind(NhlName(am_ids).eq."right_axis")  
  if(.not.ismissing(index)) then                                             
    NhlRemoveAnnotation(plot@contour,am_ids(index))                          
  end if  
  draw(plot)


;RH plot:
  resL@cnLevelSelectionMode = "AutomaticLevels"   
  plot = gsn_csm_pres_hgt (wks,slice2(:,0:ntim-1),resL) ;
                                                               
  ; to remove the "height" label on the right, we have to go through several
; steps. The label is not an axis string but extra text placed there by
; the plot template
  
  getvalues plot@contour                                                     
  "pmAnnoManagers" : am_ids                                                
  end getvalues   

  ;index = ind(NhlName(am_ids).eq."right_axis")  
  if(.not.ismissing(index)) then                                             
    NhlRemoveAnnotation(plot@contour,am_ids(index))                          
  end if  


  draw(plot)                                                                 
  frame(wks)  

 
 ;plot = gsn_csm_xy2(wks,ispan(0,30,1),slice(:,0:30),slice2(:,0:30),resL,resR) 
end
