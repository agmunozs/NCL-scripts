load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"  
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRFUserARW.ncl"

begin

ini=0
fin=8
plev=850

;****************************************************************************
;****************************************************************************
;***********FIN INTERVENCION USUARIO*****************************************


  a = addfile("wrfout.nc","r")
  lat2d = a->XLAT(0,:,:)
  lon2d = a->XLONG(0,:,:)
  grid_sz = dimsizes(lat2d)

  opt = True
  opt@returnInt = True                    
  ;locNO = wrf_user_ll_to_ij (a, minlon, maxlat, opt)
  ;locSE = wrf_user_ll_to_ij (a, maxlon, minlat, opt)


  type = "x11"
; type = "pdf"
; type = "ps"
; type = "ncgm"


  res = True
  res@tiUseMainAttributes=True
  res@tiMainFontHeightF = 0.06
  res@gsnMaximize=True
  res@Footer=False
  res@MainTitlePos="Center"
  res@InitTime=False
  ;res@gsnDraw       = False          ; do not draw
  ;res@gsnFrame      = False          ; do not advance

  wks = gsn_open_wks(type,"precip_winds")
  res@MainTitle = ""

  pltres = True
  pltres@gsnMaximize=True
  pltres@vcGlyphStyle = "CurlyVector"

  mpres = True
  mpres@mpGeophysicalLineColor = "Black"
  mpres@mpNationalLineColor    = "Black"
  mpres@mpGridLineColor        = "Black"
  mpres@mpLimbLineColor        = "Black"
  mpres@mpPerimLineColor       = "Black"
  ;mpres@mpOutlineBoundarySets = "National";Geophysical"
  mpres@mpDataSetName = "Earth..4"
  mpres@mpDataBaseVersion = "MediumRes"
  mpres@mpOutlineOn = True 
  mpres@mpOutlineBoundarySets = "AllBoundaries" 
  mpres@gsnMaximize=True 
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; First get the variables we will need        

    ;slp = wrf_user_getvar(a,"slp",-1)    ; slp
    ;  wrf_smooth_2d( slp, 3 )            ; smooth slp
  ; Get non-convective, convective 
  ; Calculate total precipitation
    rainc    = wrf_user_getvar(a,"RAINC",-1)
    rainnc   = wrf_user_getvar(a,"RAINNC",-1)
    t2a      = wrf_user_getvar(a,"T2",-1)
    Ue       = wrf_user_getvar(a,"U",-1)
    Ve       = wrf_user_getvar(a,"V",-1)
    p        = wrf_user_getvar(a,"pressure",-1)
    rain_tot = rainc + rainnc 
printVarSummary(p)
printVarSummary(Ue)
    U        = new((/1, 1, grid_sz(0), grid_sz(1)+1/), "float")
    V        = new((/1, 1, grid_sz(0), grid_sz(1)/), "float") 
printVarSummary(U)
    U(0,:,:,:)= wrf_user_intrp3d(Ue(:,:,:,:48),p,"h",plev,0.,False)
    V(0,:,:,:)= wrf_user_intrp3d(Ve,p,"h",plev,0.,False)

    rain_tot@description = "Total Precipitation"
    t2a@description = "Average Temperature 2 m"
    U@description = "Surface Winds"

    ;rain_tot@description = "PRECIPITACION TOTAL"
    ;t2a@description = "TEMPERATURA MEDIA A 2 m"
    ;U@description = "VIENTOS"
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


   ; Plotting options for Precipitation
     opts_r = res                        
     opts_r@UnitLabel            = "mm"
     opts_r@cnFillOn                 = True             ; turn on color
     opts_r@cnLinesOn                = False;           ; turn off contour lines
     opts_r@cnLevelSpacingF          = 10.              ; contour interval
     opts_r@cnLevelSelectionMode     = "ManualLevels"
     opts_r@cnMinLevelValF           = 0.
     opts_r@cnMaxLevelValF           = 120.

     mpres@mpProjection          = "Lambert"
     ;mpres@mpLimitMode           = "Corners"  
     ;mpres@mpLeftCornerLatF      = lat2d(locNO(1),locNO(0))
     ;mpres@mpLeftCornerLonF      = lon2d(locNO(1),locNO(0))
     ;mpres@mpRightCornerLatF     = lat2d(locSE(1),locSE(0))
     ;mpres@mpRightCornerLonF     = lon2d(locSE(1),locSE(0))

     opts_r@gsnMaximize          =True
     opts_r@lbOrientation        ="vertical"
     opts_r@lbTitleOn            = False
     opts_r@pmLabelBarSide       ="right"
     opts_r@pmLabelBarOrthogonalPosF= 0.05
     opts_r@cnInfoLabelOn        = False
     opts_r@cnConstFLabelOn      = False
     opts_r@cnFillOn             = True
     ;pltres@gsnLeftString        = "Precipitacion Acumulada (mm)"
     ;pltres@gsnRightString       = "CMC"
     pltres@NoHeaderFooter     = True
     ;pltres@vcRefLengthF            = 0.050            ; ref vec length
    ; pltres@vcRefAnnoString2On = True   
    pltres@vcRefMagnitudeF = 3                ; define vector ref mag
    pltres@vcRefLengthF              = 0.040           ; define length of vec ref
     pltres@vcRefAnnoOn = "True"
     ;pltres@vcRefAnnoOrthogonalPosF   = -0.75            ; move ref vector
     ;pltres@vcRefAnnoArrowLineColor   = "blue1"         ; change ref vector color
     ;pltres@vcRefAnnoArrowUseVecColor = True
     ;pltres@vcRefAnnoExplicitMagnitudeF = 2.0 
     ;pltres@vcLabelsOn=True
     ;pltres@vcRefAnnoSide = "top"
      ;pltres@vcGlyphStyle= "WindBarb"
      gsn_define_colormap(wks,"WhBlGrYeRe");
      gsn_define_colormap(wks2,"WhBlGrYeRe");


   ; Total Precipitation (color fill)
     raint = rain_tot(fin,:,:)-rain_tot(ini,:,:)
     contour_tot = wrf_contour(a,wks, raint, opts_r)
 
   ; Precipitation Tendencies 
     ;opts_r@SubFieldTitle = "from " + times(it-2) + " to " + times(it)
     opts_r@NoHeaderFooter             = True
     opts_r@cnFillOn = False
     opts_r@cnLineColor = "Red4"
     ;contour_prm = wrf_contour(a,wks,rain_tot(per,:,:),opts_r)   ; con (red lines)
     delete(opts_r)

     opts_t = res
     opts_t@NoHeaderFooter             = True
     opts_t@UnitLabel            = "~S~o~N~C"
     opts_t@cnFillOn                 = True             ; turn on color
     opts_t@cnLinesOn                = False;           ; turn off contour lines
     opts_t@cnLevelSpacingF          = 1.              ; contour interval
     opts_t@cnLevelSelectionMode     = "ManualLevels"
     opts_t@cnMinLevelValF           = 5.
     opts_t@cnMaxLevelValF           = 38.
     opts_t@lbOrientation        ="vertical"
     opts_t@lbTitleOn            = False
     opts_t@pmLabelBarSide       ="right"
     opts_t@pmLabelBarOrthogonalPosF=0.05
     opts_t@cnInfoLabelOn        = False
     opts_t@cnConstFLabelOn      = False
     opts_t@cnFillOn             = True
;uv1 = gsn_csm_vector_map(wks,U(per,:,:),V(per,:,:),pltres)
;uv = gsn_csm_vector_map(wks2,U(per,:,:),V(per,:,:),pltres)
     uv1 = wrf_vector(a,wks,U(fin,:,:),V(fin,:,:),pltres)
     t2_con = wrf_contour(a,wks2,t2a(fin,:,:),opts_t)
     uv = wrf_vector(a,wks2,U(fin,:,:),V(fin,:,:),pltres)
;************************************************
;  txres               = True                     ; text mods desired
;  txres@txFontHeightF = 0.014                     ; font smaller. default big

;  gsn_text_ndc(wks,"CMC - LUZ (www.cmc.org.ve)",0.55,0.030,txres)
;  gsn_text_ndc(wks2,"CMC - LUZ (www.cmc.org.ve)",0.55,0.030,txres)
;************************************************
   ; MAKE PLOTS                                       

     ;Total Precipitation 
     plot = wrf_map_overlays(a,wks,(/contour_tot,uv1/),pltres,mpres)
;plot = gsn_csm_contour_map_overlay(wks,contour_tot,uv1,pltres,mpres)
    ;T a 2m
    plot2 = wrf_map_overlays(a,wks2,(/t2_con,uv/),pltres,mpres)

    ;UV superficial
    ;plot = wrf_map_overlays(a,wks3,uv,pltres,mpres)


;draw(plot)
;frame(wks)
;frame(wks2)

end
