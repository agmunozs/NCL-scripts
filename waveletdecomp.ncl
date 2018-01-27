load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl" 
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl" 
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"

begin
  fname  = "PorcentajeLemnaPromedio.txt"
  data   = asciiread(fname,-1,"string")
  fecha  = stringtoint(str_get_cols(data,0,7))
  x      = stringtofloat(str_get_cols(data,8,13))
  x@_FillValue = -99999.99

  fname2 = "PorcentajeLemnaPromedioNPCI.txt"
  data2  = asciiread(fname2,-1,"string")
  x2     = stringtofloat(str_get_cols(data2,8,13))
  x2@_FillValue = -99999.99

  Ndim   = dimsizes(x)   ; number of elements 
  N      = Ndim(0)
  time0  = yyyymmdd_to_yyyyfrac(fecha, 0.5)
  time   = fspan(time0(0),time0(N-1),N)

  peri0  = (/0,15/);15
  peri1  = (/7,19/);19
  np     = dimsizes(peri0)

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



  x  = dtrend_n(x,True,0)
print(x@slope)
print(x@y_intercept)
siglvl = 0.05
  aveX = avg (x)
  aveY = avg (x@y_intercept+time*x@slope)
  varX = variance (x)
  varY = variance (x@y_intercept+time*x@slope)
  sX   = dimsizes (x)        ; X and Y can be of
  sY   = dimsizes (x@y_intercept+time*x@slope)        ; different sizes
  iflag= False               ; population variance similar
  prob = ttest(aveX,varX,sX, aveY,varY,sY, iflag, False) 
print("Si <0.05 es estad’sticamente significativo NDVI al 95%: "+prob)

  x2 = dtrend_n(x2,True,0)
print(x2@slope)
print(x2@y_intercept)
  aveX = avg (x2)
  aveY = avg (x2@y_intercept+time*x2@slope)
  varX = variance (x2)
  varY = variance (x2@y_intercept+time*x2@slope)
  sX   = dimsizes (x2)        ; X and Y can be of
  sY   = dimsizes (x2@y_intercept+time*x2@slope)        ; different sizes
  iflag= False               ; population variance similar
  prob = ttest(aveX,varX,sX, aveY,varY,sY, iflag, False) 
print("Si <0.05 es estad’sticamente significativo NDVI al 95%: "+prob)

  ;asciiwrite ("Maracaibo.dat", xma)


;************************************
; compute wavelet
;************************************

  w = wavelet(x,mother,dt,param,s0,dj,jtot,npad,noise,isigtest,siglvl,nadof)
  w2= wavelet(x2,mother,dt,param,s0,dj,jtot,npad,noise,isigtest,siglvl,nadof)

;************************************
; create coodinate arrays for plot
;************************************
  power            = onedtond(w@power,(/jtot,N/))
  power!0          = "period"                        ; Y axis
  power&period     = w@period                      ; convert period to units of years


  power!1          = "time"                          ; X axis
  power&time       = time

  power@long_name  = "Power Spectrum"
  power@units      = "1/unit-freq"

  filter  = new((/np,N/),float)
  filter2 = new((/np,N/),float)
print(w@period)


  scale   = w@scale
  scale!0 = "period"
  scale&period = power&period
printVarSummary(scale)
  Cdelta = w@cdelta


  power2            = onedtond(w2@power,(/jtot,N/))
  power2!0          = "period"                        ; Y axis
  power2&period     = w2@period                      ; convert period to units of years


  power2!1          = "time"                          ; X axis
  power2&time       = time

  power2@long_name  = "Power Spectrum"
  power2@units      = "1/unit-freq"


  scale2   = w2@scale
  scale2!0 = "period"
  scale2&period = power2&period
  Cdelta2 = w2@cdelta
 

  do k=0,np-1
    frei=peri0(k)
    fref=peri1(k)
print(frei)
print(fref)
   do it=0,N-1
    filter(k,it)  =0.0
    filter2(k,it) =0.0
    do fre=frei,fref
      filter(k,it) =filter(k,it)+power(fre,it)/sqrt(scale(fre))
      filter2(k,it)=filter2(k,it)+power2(fre,it)/sqrt(scale2(fre))
    end do
    filter(k,it)  = dj*sqrt(dt)/Cdelta/w@psi0*filter(k,it)
    filter2(k,it) = dj*sqrt(dt)/Cdelta2/w2@psi0*filter2(k,it)
   end do
  end do

;********************************************************************************
; PLOTS
;********************************************************************************
 

  wks = gsn_open_wks("png","ndvi_filters")               ; open ps file
  gsn_define_colormap(wks,"BlAqGrYeOrReVi200")    ; choose colormap

			; overlay probability plot onto power plot
 
;***********Interannual
  ress = True
  ress@xyDashPatterns = (/0,0,0,0/)
  ress@xyLineColors   = (/"red","blue","green","brown"/)
  ress@xyLineThicknesses = (/4.0,4.0,4.0,4.0/)
  ress@pmLegendDisplayMode = "Always"
  ress@pmLegendOrthogonalPosF = -1.2
  ress@pmLegendParallelPosF   =  0.82
  ress@pmLegendWidthF         =  0.22
  ress@pmLegendHeightF        =  0.1
  ress@lgLabelFontHeightF     =  0.015
  ress@lgPerimOn              = False                  ; turn off box around
  ress@gsnFrame = False
  ress@gsnDraw = False
  ress@vpHeightF           = .27
  ress@vpWidthF            = .85
  ress@tmXBPrecision       = 4
  ress@gsnYRefLine   = 0
  ress@gsnYRefLineThicknessF = 4.
  ress@tmXBTickStartF      = 2005.
  ress@tmXBTickEndF        = 2012.
  ress@tiYAxisFontHeightF  = 0.02
  ress@tmXBLabelsOn         = False              ; do not draw bottom labels
  ress@tmXBOn               = False   

plot = new(3,graphic)
;***********Interdecadal
  ress@xyExplicitLegendLabels = (/"NDVI","NPCI","Ensemble"/)
  ress@tiYAxisString = "Raw Anom. (%)"
  ress@tiXAxisString = ""
  ress@vpYF         = 0.95
  ress@vpXF         = 0.125
  ress@tmXTOn       = True
  ;ress@tmXTLabelsOn = False
  ress@tmXBOn       = False

  plot(0)= gsn_csm_xy(wks,power&time,(/x,x2,(x+x2)*0.5/),ress)
;,(x+x2)*0.5
;ress@xyLineThicknesses = (/1.5,1.5,1.0,4.0/)


  ress@vpYF         = 0.65
  ress@tmXTOn       = True
  ;ress@tmXTLabelsOn = False
  ress@tmXBOn       = True
  ress@tiYAxisString = "Weekly Anom. (%)"
  plot(1)= gsn_csm_xy(wks,power&time,(/dtrend_n(filter(0,:),False,0),dtrend_n(filter2(0,:),False,0),(dtrend_n(filter(0,:),False,0)+dtrend_n(filter2(0,:),False,0))*0.5/),ress)
;(dtrend_n(filter(0,:),False,0)+dtrend_n(filter2(0,:),False,0))*0.5


  ress@vpYF         = 0.35
  ress@tmXTOn       = False
  ress@tmXTLabelsOn = False

;ress@xyLineThicknesses = (/4.0,4.0,4.0,4.0/)
  ress@tiYAxisString = "Monthly Anom. (%)"
  ress@tmXBLabelsOn         = True              ; do not draw bottom labels
  ress@tmXBOn               = True  
plot(2)= gsn_csm_xy(wks,power&time,(/dtrend_n(filter(1,:),False,0),dtrend_n(filter2(1,:),False,0),(dtrend_n(filter(1,:),False,0)+dtrend_n(filter2(1,:),False,0))*0.5/),ress)

  pres = True
  pres@gsnMaximize = True
  pres@gsnPaperOrientation = "portrait"
  pres@gsnPanelBottom = 0.05
  pres@gsnPanelFigureStrings = (/"a","b","c"/)
  ;pres@gsnPanelBottom   = 0.2		
  ;pres@gsnPanelTop      = 0.9
  pres@gsnPanelYWhiteSpacePercent = 0.
  draw(wks)
  frame(wks)
  ;gsn_panel(wks,plot,(/3,1/),pres)
end
