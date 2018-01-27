load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"

begin
;----------------------------------------------------------------------
var  = "PR" ; PR or T
 nm   = 3  ; number of members
 peri = (/"Jan-Dec","Dec-Feb","Jun-Aug"/)
 scal = (/"Inter-annual","Decadal","Trend"/)
 model= (/"MRI-CGCM3","CESM1-CAM5","MIROC5","CNRM-CM5","HadCM3","CSIRO-Mk3-6-0","MPI-ESM-LR","MPI-ESM-MR"/)
;----------------------------------------------------------------------
;Read data. First annual:
 do mo=0,dimsizes(model)-1     ;counter for model
  do m=1,nm				   ;counter for member
   do s = 0,dimsizes(scal)-1   ;counter for timescale
   
     data = readAsciiTable(ver+"_"+model(mo)+"_"+peri(0)+"_m"+m+".nc", 6, "float", 2)

;Read seasonal data:
data = readAsciiTable("PR_CESM1-CAM5_Dec-Feb_m1.nc", 6, "float", 2)

print(data)
;----------------------------------------------------------------------
; Open a PNG file and define the color table
;----------------------------------------------------------------------
  colors  = (/"white","black","dodgerblue2","firebrick","forestgreen",\
              "dodgerblue2","darkorange2","forestgreen","deeppink2",\
              "darkorange2","purple"/)

  wks = gsn_open_wks ("png","scatter")
  gsn_define_colormap(wks,colors)

;----------------------------------------------------------------------
; Set some plot resources
;----------------------------------------------------------------------
  res                    = True

  res@gsnDraw            = False       ; Don't draw plot
  res@gsnFrame           = False       ; Don't advance frame

;
; Don't use gsnMaximize. Instead, we will control 
; exactly where plot goes.
;
  res@vpXF              = 0.2
  res@vpYF              = 0.85
  res@vpWidthF          = 0.5
  res@vpHeightF         = 0.5

;---Some calculations to get min/max for Y axis.
  miny = min(avdataP) 
  maxy = max(avdataP) 
  minx = min(avdataT) 
  maxx = max(avdataT) 

  amaxy = maxy - abs(miny)
  amaxy = amaxy / 2.
  amaxy = amaxy + 25.
  
  aminy = maxy - abs(miny)
  aminy = aminy / 2.
  aminy = aminy - 25.

;---Resources for X/Y axes.
  res@trYMinF            =  aminy               ; min value on y-axis
  res@trYMaxF            =  amaxy               ; max value on y-axis
  res@trXMinF            =  0.                  ; min value on x-axis
  res@trXMaxF            =  5.                  ; max value on x-axis

;---XY marker resources
  res@xyMarkLineMode     = "Markers"
  res@xyMarkerSizeF      = 0.03
  res@xyMarkerThicknessF = 3.
  res@xyMonoMarkerColor  = False          ; Allow different colors for markers
  res@xyMarkerColors     = colors(2:)     ; It's okay to list more than
                                          ; you need here
  res@xyMarkers          = (/4,4,4,7,7,7,3,6,6/)

;---Tickmark resources
  res@tmXBMode         = "Manual"
  res@tmYLMode         = "Manual"
  res@tmXBTickSpacingF = 1.
  res@tmYLTickSpacingF = 5.

;---Title resources
  res@tiMainString       = "30 av Annual % Change P vs Delta T~C~" + \
                           "               for the WetSouth"
  res@tiXAxisString      = "Delta T"
  res@tiYAxisString      = "% Change P"
  res@tiMainFontHeightF  = 0.03

;----------------------------------------------------------------------
; Create the scatter plot, but don't draw it yet.
;----------------------------------------------------------------------
  plot = gsn_csm_xy (wks,transpose((/avdataT,avdataT/)), \
                         transpose((/avdataP,avdataP/)),res)

;----------------------------------------------------------------------
; Draw time stamp first, then the plot.
;----------------------------------------------------------------------
  infoTimeStamp(wks,0.005,"")
  draw(plot)

;----------------------------------------------------------------------
; Draw some individual labelbars.
;----------------------------------------------------------------------
  labels1  = (/"CRCM","ECP2","HRM3","MM5I","RCM3","WRFG"/)
  colors1  = (/"dodgerblue2","purple","deeppink2","firebrick",\
               "darkorange2","forestgreen"/)
  lbres                    = True          ; labelbar only resources
  lbres@vpWidthF           = 0.1           ; labelbar width
  lbres@vpHeightF          = 0.1           ; labelbar height
  lbres@lbBoxMajorExtentF  = 0.15          ; puts space between color boxes
  lbres@lbMonoFillPattern  = True          ; Solid fill pattern
  lbres@lbLabelFontHeightF = 0.015         ; font height. default is small
  lbres@lbLabelJust        = "CenterLeft"  ; left justify labels
  lbres@lbPerimOn          = False
;
; Each labelbar has just one label. This allows you to more
; easily control where the label goes.
; 
  xpos = (/0.05,  0.20,   0.35,   0.50,  0.65,   0.8/)
  do i=0,5
     lbres@lbFillColors     = colors1(i)
     lbres@lbLabelFontColor = colors1(i)
     gsn_labelbar_ndc(wks,1,labels1(i),xpos(i),0.2,lbres)
  end do

;----------------------------------------------------------------------
; Draw some markers and text.
;----------------------------------------------------------------------
  labels2  = (/"CCSM","CGCM3","GFDL","HADCM3","NCEP"/)
  markers1 = (/     4,      7,     6,       3,     2/)
  xpos2    = (/  0.15,   0.30,  0.45,    0.60,  0.75/)
  xpos3    = (/  0.18,   0.33,  0.48,    0.63,  0.78/)

  mkres               = True         ; Marker resources
  txres               = True         ; Text resources
  txres@txFontHeightF = 0.015
  txres@txJust        = "CenterLeft"
  do i = 0,4
     mkres@gsMarkerThicknessF = 3.5
     mkres@gsMarkerSizeF      = 0.02
     mkres@gsMarkerIndex      = markers1(i)

     gsn_polymarker_ndc(wks,           xpos2(i),0.23,mkres)
     gsn_text_ndc      (wks,labels2(i),xpos3(i),0.23,txres)
  end do

  frame(wks)    ; Now advance the frame.
end