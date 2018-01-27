load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
begin
;************************************************
; open file and read in variable
;************************************************
  fname  = "./prueba2.txt"
  data   = asciiread(fname,-1,"string")
  algo   = (str_get_cols(data,6,7))
  year   = stringtointeger(str_get_cols(data,6,7))
  month  = stringtoint(str_get_cols(data,3,4))
  day    = stringtoint(str_get_cols(data,0,1))
  eco    = stringtofloat(str_get_cols(data,9,14))
  eco@_FillValue=-99.9

  TIME = yyyymmdd_time(2009, 2010, "integer")
  time = yyyymmdd_to_yyyyfrac(TIME({20100731:20101109}),0.5)
  print(time+"  "+eco)

  fname2 = "Variables.txt"
  ncol   = numAsciiCol(fname2)
  print(ncol)
  data2  = readAsciiTable(fname2, ncol, "float", 0)
;  data2  = asciiread(fname2,-1,"string")
  varclim=(/"Temperatura Media","Temperatura Max","Temperatura Min","Radiacion Solar de Superficie","Radiacion Solar","Nubosidad","Humedad Especifica","Vientos Zonales","Vientos Meridionales","Precipitacion","Precipitacion convectiva"/)
  do i=1,ncol-1
  clim   = data2(:,i)
;************************************************
; calculate cross correlations
;************************************************
  mxlag = 10                                   ; set lag
; note, the max lag should not be more than N/4
  
;  ccr = esccr(eco,clim,maxlag)                   ; calc cross cor
  x   = ispan(-mxlag,mxlag,1)                     ; define x axis
  x_Lead_y = esccr(eco,clim,mxlag)
  y_Lead_x = esccr(clim,eco,mxlag)    ; switch the order of the series

     ccr = new ( 2*mxlag+1, float)    
     ccr(0:mxlag-1) = y_Lead_x(1:mxlag:-1)  ; "negative lag", -1 reverses order
     ccr(mxlag:)    = x_Lead_y(0:mxlag)     ; "positive lag"
;************************************************
; plot the correlations
;************************************************
  wks  = gsn_open_wks("png","correlacion_col"+i)                  ; open a ps plot
  res                 = True                         ; make plot mods
  res@tiMainString    = "Exito de anidacion y "+varclim(i-1) ; title
  res@tiXAxisString   = "Fecha"                        ; x-axis label
  res@tiYAxisString   = ""
  res@trXMinF = min(time)
  res@trXMaxF = max(time)
  res@xyLineColors = (/"red"/)
  res2=res
  res2@xyLineColors = (/"blue"/)
  plot = gsn_csm_xy2(wks,time,eco,clim,res,res2)  


  res@tiMainString    = "Exito de anidacion vs "+varclim(i-1) ; title
  res@tiXAxisString   = "Sesgo (dias)"                        ; x-axis label
  res@tiYAxisString   = "Correlacion"
  res@trXMinF = -mxlag
  res@trXMaxF = mxlag

  plot = gsn_xy(wks,x,ccr,res)                       ; plot correlation
;************************************************
end do
end
