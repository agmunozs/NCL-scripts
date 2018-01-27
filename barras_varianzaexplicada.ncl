; These files are loaded by default in NCL V6.2.0 and newer
; load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
; load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"

begin
 var  = "PR" ; PR or T
 nm   = 3  ; number of members
 peri = (/"Jan-Dec","Dec-Feb","Jun-Aug"/)
 scal = (/"Trend","Decadal","Inter-annual"/)
 model= (/"MRI-CGCM3","CESM1-CAM5","MIROC5","CNRM-CM5","HadCM3","CSIRO-Mk3-6-0","MPI-ESM-LR","MPI-ESM-MR"/)

;----------------------------------------------------------------------
 varexp = new((/dimsizes(model),dimsizes(scal),nm/),"float") ; create array to hold explained variance

;Read data. First annual:
 do mo=0,dimsizes(model)-1     ;counter for model
  do m=1,nm                    ;counter for member
   do s = 0,dimsizes(scal)-1   ;counter for timescale

     data = readAsciiTable(var+"_"+model(mo)+"_"+peri(0)+"_m"+m+".tsv", 5, "float", 2)
       				
  varexp(mo,s,m-1)=dim_variance(data(:,s+2))/dim_variance(data(:,1))*100
 
 print(peri(0)+" Model:"+model(mo)+" - Timescale:"+scal(s)+" - Member:"+m)
;  print(varexp(mo,s,m-1))
  delete(data)
  end do
 end do
end do
  
;======================================================================================
; Start the graphics
;======================================================================================

  wks      = gsn_open_wks("pdf","explainedvariance_intermodel_dispersion_"+peri(0))          ; send graphics to PNG file

  sres                       = True

  sres@gsnDraw               = True
  sres@gsnFrame              = False

  sres@vpWidthF              = 0.7
  sres@vpHeightF             = 0.5
  sres@vpXF                  = .15

  sres@trXMinF               =  0.4
  sres@trXMaxF               = 8.4
  sres@trYMinF               = 0
  sres@trYMaxF               =  100

  sres@gsnXYBarChart         = True
  sres@gsnXYBarChartBarWidth = 0.15     ; change bar widths
  sres@gsnYRefLineColor      = "transparent"

  sres@tiYAxisSide		     = "Right"
  sres@tmYRLabelsOn			 = True
  sres@tmYLLabelsOn			 = False
  ;sres@tmXTLabelsOn			 = True
  sres@tmXBMode              = "Explicit"      ; explicit labels
  sres@tmXBValues            = (/1,2,3,4,5,6,7,8/)
  sres@tmXBLabels            = model
  sres@tmXBLabelAngleF 		 = 90
  sres@tmYLLabelAngleF 		 = 90
  sres@tmYRLabelAngleF 		 = 90
  sres@tmXBLabelFontHeightF  = 0.0125
  sres@tmXTLabelFontHeightF  = 0.0125
  sres@tmYLLabelFontHeightF  = 0.0125

  sres@tiMainString          = ""
  sres@tiMainFont            = "helvetica"
  sres@tiMainFontHeightF     = 0.0125
  sres@tiYAxisString         = "Explained Variance (%)"
  sres@tiYAxisFontHeightF     = 0.0125
  sres@tmYMajorGrid			  = True
  sres@tmYMajorGridThicknessF = 1.0           ; 2.0 is default
  sres@tmYMajorGridLineDashPattern = 2  
  sres@tmYMinorGrid			  = True
  sres@tmYMinorGridThicknessF = 1.0           ; 2.0 is default
  sres@tmYMinorGridLineDashPattern = 2  
  sres@tmYLOn					= False
  sres@tmXTOn					= False
  
  sres@gsnRightString        = ""

    
  plot1 = new((/2,dimsizes(model)/),graphic)
  plot2 = plot1
  plot3 = plot1
  do gg = 0,dimsizes(model)-1     ;counter for model
     sres@gsnXYBarChartColors = "red"
     sres@gsnYRefLine         = dim_avg(varexp(gg,0,:))    
     plot1(0,gg) = gsn_csm_xy(wks,.775+gg,max(varexp(gg,0,:)),sres)    
     plot1(1,gg) = gsn_csm_xy(wks,.775+gg,min(varexp(gg,0,:)),sres)    
       
     sres@gsnXYBarChartColors = (/"green"/)
     sres@gsnYRefLine         = dim_avg(varexp(gg,1,:))     

     plot2(0,gg) = gsn_csm_xy(wks,.925+gg,max(varexp(gg,1,:)),sres)
     plot2(1,gg) = gsn_csm_xy(wks,.925+gg,min(varexp(gg,1,:)),sres)   

     sres@gsnXYBarChartColors = (/"blue"/)
     sres@gsnYRefLine         = dim_avg(varexp(gg,2,:))     
     plot3(0,gg) = gsn_csm_xy(wks,1.075+gg,max(varexp(gg,2,:)),sres)
     plot3(1,gg) = gsn_csm_xy(wks,1.075+gg,min(varexp(gg,2,:)),sres)   
      
  end do
  
;---Create labelbars
  lbres                    = True          ; labelbar only resources
  lbres@vpWidthF           = 0.09           ; labelbar width
  lbres@vpHeightF          = 0.5           ; labelbar height
  lbres@lbBoxMajorExtentF  = 0.36          ; puts space between color boxes
  lbres@lbMonoFillPattern  = True          ; Solid fill pattern
  lbres@lbFillColors       = (/"red","green","blue"/)
  lbres@lbLabelFontHeightF = 0.0125         ; font height. default is small
  ;lbres@lbLabelAlignment   = "ExternalEdges"
  lbres@lbLabelJust        = "BottomLeft"  ; left justify labels
  lbres@lbLabelAngleF      = 90
  lbres@lbPerimOn          = False
  lbres@lgPerimColor       = "white"
  lbres@gsnOrientation 	   = "landscape"
  lbres@gsnMaximize		   = True

  labels = (/"","",""/)
  gsn_labelbar_ndc(wks,3,labels,0.92,0.8,lbres)    ; draw right labelbar column
	
;************************************************
; Attach text to plot using plot coordinates.
;************************************************
  txres               = True                     ; text mods desired
  txres@txFontHeightF = 0.0125                   ; font smaller. default big
  txres@txAngleF	  = 90
  gsn_text_ndc(wks,"Trend",0.99,0.385,txres) 
  gsn_text_ndc(wks,"Decadal",0.99,0.55,txres)  
  gsn_text_ndc(wks,"Inter-annual",0.99,0.715,txres) 

  frame(wks)
  delete(wks)
  
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
;Now seasons 
;IMPORTANT: cross-year season files have a different number of columns!!--------

;Read data. 
do se=1,dimsizes(peri)-1	   ;two seasons; peri(0) is annual, already plotted.
 do mo=0,dimsizes(model)-1     ;counter for model
  do m=1,nm                    ;counter for member
   do s = 0,dimsizes(scal)-1   ;counter for timescale
   
 if(se.eq.1)
  data = readAsciiTable(var+"_"+model(mo)+"_"+peri(se)+"_m"+m+".tsv", 6, "float", 2)
  varexp(mo,s,m-1)=dim_variance(data(:,s+3))/dim_variance(data(:,2))*100
  else
  data = readAsciiTable(var+"_"+model(mo)+"_"+peri(se)+"_m"+m+".tsv", 5, "float", 2)
  varexp(mo,s,m-1)=dim_variance(data(:,s+2))/dim_variance(data(:,1))*100
 end if
 
 print(peri(se)+" Model:"+model(mo)+" - Timescale:"+scal(s)+" - Member:"+m)
;  print(varexp(mo,s,m-1))
  delete(data)
  end do
 end do
end do
  
;======================================================================================
; Start the graphics
;======================================================================================

  wks      = gsn_open_wks("pdf","explainedvariance_intermodel_dispersion_"+peri(se))          ; send graphics to PNG file

  sres                       = True

  sres@gsnDraw               = True
  sres@gsnFrame              = False

  sres@vpWidthF              = 0.7
  sres@vpHeightF             = 0.5
  sres@vpXF                  = .15

  sres@trXMinF               =  0.4
  sres@trXMaxF               = 8.4
  sres@trYMinF               = 0
  sres@trYMaxF               =  100

  sres@gsnXYBarChart         = True
  sres@gsnXYBarChartBarWidth = 0.15     ; change bar widths
  sres@gsnYRefLineColor      = "transparent"

  sres@tiYAxisSide		     = "Right"
  sres@tmYRLabelsOn			 = True
  sres@tmYLLabelsOn			 = False
  ;sres@tmXTLabelsOn			 = True
  sres@tmXBMode              = "Explicit"      ; explicit labels
  sres@tmXBValues            = (/1,2,3,4,5,6,7,8/)
  sres@tmXBLabels            = model
  sres@tmXBLabelAngleF 		 = 90
  sres@tmYLLabelAngleF 		 = 90
  sres@tmYRLabelAngleF 		 = 90
  sres@tmXBLabelFontHeightF  = 0.0125
  sres@tmXTLabelFontHeightF  = 0.0125
  sres@tmYLLabelFontHeightF  = 0.0125

  sres@tiMainString          = ""
  sres@tiMainFont            = "helvetica"
  sres@tiMainFontHeightF     = 0.0125
  sres@tiYAxisString         = "Explained Variance (%)"
  sres@tiYAxisFontHeightF     = 0.0125
  sres@tmYMajorGrid			  = True
  sres@tmYMajorGridThicknessF = 1.0           ; 2.0 is default
  sres@tmYMajorGridLineDashPattern = 2  
  sres@tmYMinorGrid			  = True
  sres@tmYMinorGridThicknessF = 1.0           ; 2.0 is default
  sres@tmYMinorGridLineDashPattern = 2  
  sres@tmYLOn					= False
  sres@tmXTOn					= False
  
  sres@gsnRightString        = ""

    
  plot1 = new((/2,dimsizes(model)/),graphic)
  plot2 = plot1
  plot3 = plot1
  do gg = 0,dimsizes(model)-1     ;counter for model
     sres@gsnXYBarChartColors = "red"
     sres@gsnYRefLine         = dim_avg(varexp(gg,0,:))    
     plot1(0,gg) = gsn_csm_xy(wks,.775+gg,max(varexp(gg,0,:)),sres)    
     plot1(1,gg) = gsn_csm_xy(wks,.775+gg,min(varexp(gg,0,:)),sres)    
       
     sres@gsnXYBarChartColors = (/"green"/)
     sres@gsnYRefLine         = dim_avg(varexp(gg,1,:))     

     plot2(0,gg) = gsn_csm_xy(wks,.925+gg,max(varexp(gg,1,:)),sres)
     plot2(1,gg) = gsn_csm_xy(wks,.925+gg,min(varexp(gg,1,:)),sres)   

     sres@gsnXYBarChartColors = (/"blue"/)
     sres@gsnYRefLine         = dim_avg(varexp(gg,2,:))     
     plot3(0,gg) = gsn_csm_xy(wks,1.075+gg,max(varexp(gg,2,:)),sres)
     plot3(1,gg) = gsn_csm_xy(wks,1.075+gg,min(varexp(gg,2,:)),sres)   
      
  end do
  
;---Create labelbars
  lbres                    = True          ; labelbar only resources
  lbres@vpWidthF           = 0.09           ; labelbar width
  lbres@vpHeightF          = 0.5           ; labelbar height
  lbres@lbBoxMajorExtentF  = 0.36          ; puts space between color boxes
  lbres@lbMonoFillPattern  = True          ; Solid fill pattern
  lbres@lbFillColors       = (/"red","green","blue"/)
  lbres@lbLabelFontHeightF = 0.0125         ; font height. default is small
  ;lbres@lbLabelAlignment   = "ExternalEdges"
  lbres@lbLabelJust        = "BottomLeft"  ; left justify labels
  lbres@lbLabelAngleF      = 90
  lbres@lbPerimOn          = False
  lbres@lgPerimColor       = "white"
  lbres@gsnOrientation 	   = "landscape"
  lbres@gsnMaximize		   = True

  labels = (/"","",""/)
  gsn_labelbar_ndc(wks,3,labels,0.92,0.8,lbres)    ; draw right labelbar column
	
;************************************************
; Attach text to plot using plot coordinates.
;************************************************
  txres               = True                     ; text mods desired
  txres@txFontHeightF = 0.0125                   ; font smaller. default big
  txres@txAngleF	  = 90
  gsn_text_ndc(wks,"Trend",0.99,0.385,txres) 
  gsn_text_ndc(wks,"Decadal",0.99,0.55,txres)  
  gsn_text_ndc(wks,"Inter-annual",0.99,0.715,txres) 

  frame(wks)
  delete(wks)
end do
end