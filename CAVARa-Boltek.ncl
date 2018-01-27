;====================HACER RADAR==============;
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
begin

;************************************************
;Procedimiento modular para labels
;************************************************
procedure add_lc_labels(wks,map,minlat,maxlat,minlon,maxlon,ytickvals,yticknames)
local lat_values, nlat, lat1_ndc, lat2_ndc, lon1_ndc, lon2_ndc,slope,txres, \
lon_values, PI, RAD_TO_DEG, dum_lft, dum_rgt, dum_bot
begin
;---Pick some "nice" values for the latitude labels.
  lat_values = ytickvals
  nlat       = dimsizes(lat_values)

;---Set some text resources
  txres               = True
  txres@txFontHeightF = 0.015
  txres@txPosXF       = 0.1
  ;txres@txFontColor	      = "yellow"
;
; Loop through lat values, and attach labels to the left and
; right edges of the masked LC plot.  The labels will be 
; rotated to fit the line better.
;
  dum_lft = new(nlat,graphic)    ; Dummy array to hold attached strings.
  do n=0,nlat-1
; Add extra white space to labels.
    lat_label_lft = yticknames(n)  

;---Left label
;    txres@txAngleF = 0.0 
    dum_lft(n) = gsn_add_text(wks,map,lat_label_lft,minlon,ytickvals(n),txres)

  end do

end

;CODIGO PRINCIPAL

;************************************************
; Apertura y lectura
;************************************************
  filename = asciiread("lectura",-1,"string")
  archi  = filename(0)+".txt"
  ;archi  = "20110913".txt"
  ncol   = numAsciiCol(archi)
  ;data   = readAsciiTable(archi, ncol, "float", "No_FillValue")
  ;data@_FillValue = 9.96921e+36
  x2 = asciiread(archi,-1,"float")
  data = onedtond(x2(0:), (/numAsciiRow(archi),5/))
  ;print (data)

  tipo   = data(:,3)
  polari = data(:,4) 
  dsizes = dimsizes(tipo)
  angul  = data(:,1)
  radios = data(:,2)
  max_rad = max(radios)
  angles = angul-360     ; fix to make angles monotonic

  tip    = (/"TN","IN"/)
  tipcol = (/"red","blue"/)
  pol    = (/"+","-"/)
;CARAGAR PNG
 wks_type = "png"
 wks_type@wkWidth = 1150
 wks_type@wkHeight = 1150

  wks = gsn_open_wks(wks_type,"radar")
  gsn_define_colormap(wks,"gui_default")

; rescale amplitudes so they are in the range 0-90N
; method: work from max plotted amplitude from yticks, divide by 90
; then divide all amplitudes by same factor, then do 90- answer to make outside ring 90N
yticknames=(/0,100,200,300,400,500/)
max_amp=max(yticknames)
fact=(max_amp/90.0)
radiosf=90-(radios/fact)
ytickvals=yticknames/fact	; these are the ytick values
yticknames=yticknames(::-1)
;---Set some resources for a polar map.
  res              = True

  res@gsnMaximize  = True                 ; Maximize plot in frame
  res@gsnDraw      = False                ; Don't draw plot 
  res@gsnFrame     = False                ; Don't advance frame (page)

  res@gsnTickMarksOn  = True
  res@tiMainString  = "Descargas El"+"e:H-13V2F35:B:FV-2H3:"+"ctricas "+filename(0)

  res@gsnPolar     = "NH"                 ; North hemisphere

  res@mpFillOn     = False                ; Turn off map fill
  res@mpOutlineOn  = False                ; Turn off map outlines

  res@mpCenterLonF = -90                  ; Rotate the map
  res@mpCenterLatF = 90                  ; Rotate the map
  res@mpGridLatSpacingF = 90/5.
;---Create a polar map with the given resources.
  map = gsn_csm_map_polar(wks,res)

;---Create a resource list for the markers.
  mkres = True

;---Arrays for attaching two sets of markers
  dum_hollow = new(dsizes,graphic)

  mkres@gsMarkerThicknessF = 3.0      ; Twice as thick
    
    txres               = True
    txres@gsMarkerSizeF = 0.02
    txres@txPosXF       = 0.1

    menos = NhlNewMarker(wks, "-", 18, 0.0, 0.075, 1.3125, 0.95, 0.0)
  do i=0,dsizes-1

    if (tipo(i).eq."0") then
     txres@gsMarkerColor      = tipcol(0)
     if (polari(i).eq."0") then
       txres@gsMarkerIndex =2
      else
       txres@gsMarkerIndex =menos
     end if
     else
     txres@gsMarkerColor      = tipcol(1)
     if (polari(i).eq."0") then
       txres@gsMarkerIndex =2
      else
       txres@gsMarkerIndex =menos
     end if
    end if
    ;dum_hollow(i) = gsn_add_text(wks,map,tipof,angul(i),radiosf(i),txres)
    gsn_polymarker(wks,map,angul(i),radiosf(i),txres)
  end do
;
; add y-axis labels
minlat=0
maxlat=90
minlon=120
maxlon=105
add_lc_labels(wks,map,minlat,maxlat,minlon,maxlon,ytickvals,yticknames)
;==============================================
; create xy plot
;=============================================
; x   = dim_avg(chi)                    ; average chi
; x!0 = "chi"                           ; remove warning message
; y   = ispan(0,dimsizes(chi&time)-1,1) ; create y-axis
; y!0 = "time"   

; xyres                 = True         ; xy plot mods desired
; xyres@vpWidthF        = .20          ; set width of second plot
; xyres@tmXBMinorOn     = False        ; no minor tickmarks
; xyres@tmXBLabelStride = 2            ; label stride

; xyres@gsnDraw         = False        ; don't draw yet
; xyres@gsnFrame        = False        ; don't advance frame yet

; xyres@gsnCenterString = "Zonal Ave"  ; add title
; xyres@txFontHeightF   = .015         ; change font height

; plot = gsn_csm_xy(wks, x,y,xyres)
;==============================================
;  attach plots
;=============================================
; newplot = gsn_attach_plots(base_plot,(/plot/),hres,xyres)
; draw(base_plot)
; frame(wks)
;---Drawing the map plot will draw all the attached markers too.
  draw(map)
  frame(wks)
end

 ;system("composite -dissolve 20 -gravity Center zoom.png radar.png Boltek.png")
 ;system("convert -trim Boltek.png Boltek.png")
