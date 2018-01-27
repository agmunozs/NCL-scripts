;*****************************************************
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"  
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"  
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"  
;*****************************************************
begin
 var  = "TS" ; PR or TS
 varname = "Temperature"
 nm   = 3  ; number of members
 sd	  = 2. ; how many standard deviations to mask (uncertainty threshold)
 peri = (/"Jan-Dec","Dec-Feb","Jun-Aug"/)
 scal = (/"Inter-annual","Decadal","Trend"/)
 model= (/"MRI-CGCM3","CESM1-CAM5","MIROC5","CNRM-CM5","HadCM3","CSIRO-Mk3-6-0","MPI-ESM-LR","MPI-ESM-MR"/)

 wks = gsn_open_wks("pdf","models_tdecomp"+varname)            ; open png file
 gsn_define_colormap(wks,"gui_default")        				; choose colormap
 plot = new(9, graphic)                        				; create graphical array

 k=0  						;counter for plots
 
 do mo=0,dimsizes(model)-1     ;counter for model
  ;do m=1,nm				   ;counter for member
   do s = 0,dimsizes(scal)-1   ;counter for timescale
    do p=0,dimsizes(peri)-1    ;counter for period

 print("Working on: Model:"+model(mo)+" Scale:"+scal(s)+" Period:"+peri(p))
 
  diri = "/Users/agms/Documents/Angel/Tesistas/Marling_Master/"   ; input directory
  ;fils = systemfunc ("ls "+diri+var+"_"+model(mo)+"_"+scal(s)+"_"+peri(p)+"_m"+m+".nc") ; files
  fils = systemfunc ("ls "+diri+var+"_"+model(mo)+"_"+scal(s)+"_"+peri(p)+"*.nc") ; files


  f    = addfiles (fils, "r")   

  ListSetType (f, "join")       ; join
  ratio    = f[:]->ratio        ; read field from all files  
  
  mmm = dim_avg_n_Wrap( ratio, 0 )   ; average on all members
  mmd = dim_stddev_n_Wrap( ratio, 0 )   ; STDV on all members
  mmm = mask(mmm, mmd.ge.sd,False)      ;mask all values >= SD threshold
  
  if (any(isnan_ieee(mmm))) then ;manage missing values
      value = 1.e20
      replace_ieeenan (mmm, value, 0)
      mmm@_FillValue = value
  end if
  
  ;printVarSummary (mmm)
  mmm@long_name = peri(p)
  
;***************************************
; create individual plots
;***************************************  

  res                      = True               ; plot mods desired

  res@cnFillOn             = True               ; color contours
  res@cnLinesOn            = False              ; default is True
  res@cnLineLabelsOn       = False              ; default is True
  res@cnFillMode           = "RasterFill"    ; Raster Mode

;---This resource not needed in NCL V6.1.0
  res@gsnSpreadColors      = True               ; use full range of colormap
  res@gsnStringFontHeightF = 0.022
  res@tmXBLabelFontHeightF = 0.018
  res@tmYLLabelFontHeightF = 0.018
 

  res@gsnDraw              = False              ; don't draw yet
  res@gsnFrame             = False              ; don't advance frame yet
  res@lbLabelBarOn         = False           ; turn off individual lb's
  ;res@mpProjection	   = "Mercator"
  res@gsnAddCyclic	   = False
  res@mpLimitMode 	   = "LatLon"

  res@mpMinLatF        = 4.25
  res@mpMaxLatF        = 32.75
  res@mpMinLonF        = -117.75
  res@mpMaxLonF        = -52.25

  res@cnLevelSelectionMode = "ManualLevels"     ; set manual contour levels
  res@cnMinLevelValF       = 0.                ; set min contour level
  res@cnMaxLevelValF       = 100.0                ; set max contour level
  res@cnLevelSpacingF      = .5              ; set contour spacing
  res@tmXBTickSpacingF     = 10.
  res@cnMissingValFillColor = "white"
  res@cnMissingValFillPattern = "SolidFill"
  res@gsnPaperOrientation = "landscape"


  res@gsnCenterString      = scal(s)
  plot(k)  = gsn_csm_contour_map_ce(wks,mmm,res)
  k=k+1
  ;end do
  end do        ;period
  delete(f)
  delete(ratio)
 end do			;timescale
;************************************************
; create panel plot
;************************************************

  resP                  = True                   ; modify the panel plot
  resP@txString         = "Explained Variance - "+model(mo)+" : "+varname   ; plot title
  resP@gsnMaximize      = True                   ; make ps, eps, pdf, ... large
  resP@gsnPanelLabelBar = True                   ; add common label bar
  resP@lbOrientation   = "Vertical"                ; vertical label bar
  resP@lbLabelFontHeightF = 0.01
  resP@lbLabelStride	  = 10
  resP@lbBoxSeparatorLinesOn= False
  resP@gsnPaperOrientation = "landscape"
  gsn_panel(wks,plot,(/3,3/),resP)               ; now draw as one plot

k=0  						;counter for plots
end do			; model

end
