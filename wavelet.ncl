load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl" 
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl" 
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"

begin
  initm  = 0
  fname  = "PorcentajeLemnaPromedio.txt"
  data   = asciiread(fname,-1,"string")
  fecha  = stringtoint(str_get_cols(data,0,7))
  x      = stringtofloat(str_get_cols(data,8,13))
  x@_FillValue = -99999.99

  Ndim   = dimsizes(x)   ; number of elements 
  N      = Ndim(0)
  time0  = yyyymmdd_to_yyyyfrac(fecha, 0.5)
  time   = fspan(time0(0),time0(N-1),N)

print("Fechas inicial y final: "+fecha(0)+" - "+fecha(N-1))
print("Tiempos inicial y final:"+time(0)+ " - "+time(N-1))

  mother  = 0
  param   = 6.0
  dt      = 2./365.    ;timestep
  s0      = dt
  dj      = 0.25
  jtot    = 1+floattointeger(((log10(N*dt/s0))/dj)/log10(2.))
  npad    = N
  nadof   = 0
  noise   = 1
  siglvl  = .05
  isigtest= 0

;printVarSummary(x)



 if(.not.all(ismissing(x)))


  x  = dtrend_n(x,False,0)
  ;asciiwrite ("Maracaibo.dat", xma)


;************************************
; compute wavelet
;************************************

  w= wavelet(x,mother,dt,param,s0,dj,jtot,npad,noise,isigtest,siglvl,nadof)

;************************************
; create coodinate arrays for plot
;************************************
  power            = onedtond(w@power,(/jtot,N/))
  power!0          = "period"                        ; Y axis
  power&period     = w@period                      ; convert period to units of years


  power!1          = "time"                          ; X axis
  power&time       = time

  power@long_name  = ""
  power@units      = ""

; compute significance ( >= 1 is significant)
  SIG              = power                            ; transfer meta data
  SIG              = power/conform (power,w@signif,0)
  SIG@long_name    = "Significance"
  SIG@units        = " "

;printVarSummary(SIG)
;*************************************
;Maxima
;*************************************  
  gws = w@gws
  gmaxi=(local_max_1d(gws, False,  0.25, 0)) ;for now only 1 max
  gdim=dimsizes(gmaxi)
  ;print(gmaxi)
  gmax=0.
  do j=0,gdim-1
    indice=ind(gws.eq.gmaxi(j))
    ;We test for values COI-acceptable and find absolute maximum 
    if(w@period(indice).le.600. .and. w@period(indice).ge.gmax)   ;here 180 is my max period due to the COI
      gmax=w@period(indice)
      mxindi=indice
    end if
  end do
  delete(gmaxi)


    ;We now test for statistical significance:
    kont=0
    do it=0,N-1
     if (SIG(mxindi,it).ge.1.)
      kont=kont+1
     end if
    end do
    ;print(gmax+" "+SIG(mxindi,:)+" "+kont)
    ;If >=1% is statistically significant, then define maxscale for that cell.
    ;if (int2flt(kont)/int2flt(N).ge.0.)
    ;  maxscale(ila,ilo)=gmax;/12 ;transform to years
    ;end if

   end if





;********************************************************************************
; PLOTS
;********************************************************************************
 

  wks = gsn_open_wks("png","ndvi_wavelet")               ; open ps file
  gsn_define_colormap(wks,"BlAqGrYeOrReVi200")    ; choose colormap

  res                     = True                  ; plot mods desired
  res@gsnDraw             = False                 ; Do not draw plot
  res@gsnFrame            = False                 ; Do not advance frome
  res@cnFillOn            = True                  ; turn on color
  ;res@cnFillMode          = "RasterFill"          ; turn on raster mode
  ;res@cnRasterSmoothingOn = True                  ; turn on raster smoothing
  res@pmLabelBarOrthogonalPosF = 0.83
  res@cnLinesOn           = False                 ; turn off contour lines
  res@cnLineLabelsOn      = False
  res@cnInfoLabelOn       = False
  res@gsnSpreadColors     = True                  ; use full colormap
  res@trYReverse          = True                  ; reverse y-axis
  res@tmYLMode = "Explicit"
  res@tmYLValues = (/.005,0.038,0.16,0.33,0.66,1.,2.,4./)
  res@tmYLLabels = (/".005","0.038","0.16","0.33","0.66","1.","2.","4."/)
  res@tmLabelAutoStride   = True
  res@vpHeightF           = .5                    ;
  res@vpWidthF            = .7
  ;res@cnLevelSelectionMode = "ExplicitLevels"       ; set manual contour levels
  ;res@cnLevels = (/0.5,1.,2.,4./)
  res@gsnStringFontHeightF = 0.023
  res@tiYAxisFontHeightF   = res@gsnStringFontHeightF
  res@gsnCenterString       = "NDVI Wavelet Power Spectrum"
  ;res@gsnLeftString       = "NDVI: Lake Maracaibo"
  res@tiYAxisString       = "Period (yr)"
  res@cnLevelSelectionMode = "ManualLevels" ; set manual contour levels
  res@cnMinLevelValF      = 10.00        ; set min contour level
  res@cnMaxLevelValF      = 450.00        ; set max contour level
  res@cnLevelSpacingF     = 10.00        ; set contour spacing
  res@lbLabelStride       = 4
  res@lbOrientation        = "vertical"
  res@tmXBTickStartF      = 2005
  res@tmXBTickEndF        = 2012
  res@tmXBPrecision       = 4

  res2 = True                            ; res2 probability plots
  res2@trYReverse          = True
  res2@tmYLMode = "Explicit"
  res2@tmYLValues = (/.005,0.038,0.16,0.33,0.66,1.,2.,4./)
  res2@tmYLLabels = (/".005","0.038","0.16","0.33","0.66","1.","2.","4."/)
  res2@gsnDraw             = False       ; Do not draw plot
  res2@gsnFrame            = False       ; Do not advance frome
  res2@cnLevelSelectionMode = "AutomaticLevels" ; set manual contour levels
  ;res2@cnMinLevelValF      = 0.00        ; set min contour level
  ;res2@cnMaxLevelValF      = 2.00        ; set max contour level
  ;res2@cnLevelSpacingF     = 1.00        ; set contour spacing
  res2@cnInfoLabelOn       = False
  res2@cnLinesOn           = True     ; do not draw contour lines
  res2@cnLineLabelsOn      = False       ; do not draw contour labels
  res2@cnFillScaleF        = 0.5         ; add extra density
  res2@gsnLeftString = ""
  res2@gsnRightString = ""
  res2@tmXBTickStartF      = 2005
  res2@tmXBTickEndF        = 2012
  res2@tmXBPrecision       = 4

  plot = new(3,graphic)
  plot(1) = gsn_csm_contour(wks,power,res)
  plot(1) = ShadeCOI(wks,plot(1),w,power&time,False)
  iplot = gsn_csm_contour(wks,SIG,res2)
  iplot = ShadeGtContour(iplot,1.,17)
  overlay(plot(1),iplot)				; overlay probability plot onto power plot

  scale = w@scale
  Cdelta = w@cdelta
  powernorm = power
  powernorm = power/conform(power,scale,0)
  scaleavg = dj*dt/Cdelta*dim_sum_Wrap(powernorm(time|:,{period|4./365.:13./365.}))
  ;AGMS: in the above line, we control de integrated period here: period|2.:8.

  resl = True
  resl@gsnFrame = False
  resl@gsnDraw = False
  resl@trYAxisType = "LogAxis"
  resl@trYReverse          = True                  ; reverse y-axis
  resl@tmYLMode = "Explicit"
  resl@tmYLValues = (/.005,0.038,0.16,0.33,0.66,1.,2.,4./)
  resl@tmYLLabels = (/".005","0.038","0.16","0.33","0.66","1.","2.","4."/)
  resl@gsnStringFontHeightF = res@gsnStringFontHeightF 
  resl@gsnCenterString     = "Mean Power"
  resl@xyLineThicknesses = (/3.0/)
  resl@tmYRLabelsOn  = True

  plotg = gsn_csm_xy(wks,gws,power&period,resl)

  plotc = gsn_attach_plots(plot(1),plotg,res,resl)

;***********Interannual
  ress = True
  ress@xyDashPatterns = (/0,0,0,0/)
  ress@xyLineColors   = (/"red","blue","green","yellow"/)
  ress@xyLineThicknesses = (/3.0,3.0,3.0,3.0/)
  ress@pmLegendDisplayMode = "Always"
  ress@pmLegendOrthogonalPosF = -1.1
  ress@pmLegendParallelPosF   =  0.85
  ress@pmLegendWidthF           =  0.25
  ress@pmLegendHeightF          =  0.1
  ress@lgPerimOn              = False                  ; turn off box around
  ress@gsnFrame = False
  ress@gsnDraw = False
  ress@vpHeightF           = .3
  ress@vpWidthF            = .87
  ress@tmXBPrecision       = 4
  ress@tmXBTickStartF      = 2005
  ress@tmXBTickEndF        = 2012

;***********Interdecadal
  ress@xyExplicitLegendLabels = (/"1-2 wk","0.9-2. mo","0.9-1.2 yr","1.9-2.2 yr"/)
  ress@tiYAxisString = "~F33~s~S~2"
  ress@tiXAxisString = ""
  scaleavg2 = dj*dt/Cdelta*dim_sum_Wrap(powernorm(time|:,{period|.9:1.2}))
  scaleavg3 = dj*dt/Cdelta*dim_sum_Wrap(powernorm(time|:,{period|1.9:2.2}))
  scaleavg4 = dj*dt/Cdelta*dim_sum_Wrap(powernorm(time|:,{period|0.076:0.152}))
  plot(2) = gsn_csm_xy(wks,power&time,(/scaleavg,scaleavg4,scaleavg2,scaleavg3/),ress)

;**************Serie orig
  res3   = True
  res3@gsnYRefLine   = 0
  res3@gsnAboveYRefLineColor = "Green"
  res3@gsnBelowYRefLineColor = "Blue"
  res3@gsnFrame = False
  res3@gsnDraw = False
  res3@vpHeightF           = .3
  res3@vpWidthF            = .87
  res3@tiYAxisString = "Anomaly (%)"
  res3@tiXAxisString = ""
  res3@tmXBTickStartF      = 2005
  res3@tmXBTickEndF        = 2012
  res3@tmXBPrecision       = 4

  plot(0) = gsn_csm_xy(wks,power&time,x,res3)

  pres = True
  pres@gsnMaximize = True
  ;pres@gsnPaperOrientation = "portrait"
  gsn_panel(wks,plot,(/3,1/),pres)
end
