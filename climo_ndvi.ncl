load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"  
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"  
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl" 

begin
;-------------------------------------------------------------
; User specifications
;-------------------------------------------------------------
   NC      = False                            ; create netCDF?       
   PLOT    = True                             ; sample plots?

   ymdStrt = 20050101                         ; start yyyymmdd
   ymdLast = 20091123                         ; last  
   yrStrt  = ymdStrt/10000
   yrLast  = ymdLast/10000
  
   nhar    = 4                                ; number of fourier comp      

   var     = "dvi"                            ; name of file                   
   vName   = "NDVI"                            ; name for plots

   diri    = "../datos_in/"         ; input directory   

   if (NC) then
       diro= "./"         ; output dir
       filo= var+".day.anomalies."+yrStrt+"-"+yrLast+".nc"  ; output file  
   end if

   if (PLOT) then
       wksType = "ps"
       wksName = "climo"                      ; var+"."+yrStrt+"_"+yrLast
   end if

;***********************************************************
; Read user specified time and create required yyyyddd                    
;***********************************************************
;  fname  = "../datos_in/PorcentajeLemnaPromedio.txt"
;  data   = asciiread(fname,-1,"string")
;  year   = stringtointeger(str_get_cols(data,0,3))
;  month  = stringtoint(str_get_cols(data,4,5))
;  day    = stringtoint(str_get_cols(data,6,7))
;  ndvi   = stringtofloat(str_get_cols(data,9,18))

  fname  = "../datos_out/baja.dat"
  data   = asciiread(fname,-1,"string")
  year   = stringtointeger(str_get_cols(data,12,15))
  month  = stringtoint(str_get_cols(data,16,17))
  day    = stringtoint(str_get_cols(data,18,19))
  ndvi   = stringtofloat(str_get_cols(data,25,31));9,18))

  x      = new((/dimsizes(ndvi),1,1/),typeof(ndvi))
  lo=ispan(0,0,1);dimsizes(ndvi)-1,dimsizes(ndvi))
  x(:,0,0)      =ndvi
 x!0 = "time"
 x!2 = "lon"
 x!1 = "lat"
;x&time=858
 x&lat=1
 x&lon=1
;   f       = addfile (diri+fili , "r")                          

;   time    = f->time                          ; all times on file
;   ymd     = ut_calendar(time, -2)            ; yyyymmdd
;   iStrt   = ind(ymd.eq.ymdStrt)              ; index start
;   iLast   = ind(ymd.eq.ymdLast)              ; index last 
;   delete(time)
;   delete(ymd)

;***********************************************************
; Read user specified time and create required yyyyddd                    
;***********************************************************
;   time    = f->time(iStrt:iLast)             ; time:units = "hours since"
;   TIME    = ut_calendar(time, 0)             ; type float 
;   year    = floattointeger( TIME(:,0) )
;   month   = floattointeger( TIME(:,1) )
;   day     = floattointeger( TIME(:,2) ) 
   ddd     = day_of_year(year, month, day) 
   yyyyddd = year*1000 + ddd                  ; needed for input

;***********************************************************
; Read data: short2flt                                     
;*********************************************************** 
;   x       =  short2flt( f->$var$(iStrt:iLast,:,:) )    ; convert to float 
   printVarSummary( x ) 
;***********************************************************
; Compute daily climatology: raw and then 'smoothed'  
;***********************************************************
;   xClmDay=dim_avg(x)
   xClmDay = clmDayTLL(x, yyyyddd)     ; daily climatology at each grid point
   printVarSummary(xClmDay)            
 xClmDay!0 = "time"
 xClmDay!2 = "lon"
 xClmDay!1 = "lat"
 xClmDay&lat=1
 xClmDay&lon=1

;***********************************************************
; Compute smoothed daily climatology using 'nhar' harmonics
;***********************************************************
   xClmDay_sm = smthClmDayTLL(xClmDay, nhar)
   printVarSummary(xClmDay_sm)

;***********************************************************
; Compute daily anomalies using raw and smoothed climatologies
;***********************************************************
    xAnom      = calcDayAnomTLL (x, yyyyddd, xClmDay)     
 xAnom!0 = "time"
 xAnom!2 = "lon"
 xAnom!1 = "lat"
 xAnom&time=ispan(0,dimsizes(ndvi)-1,1)
    printVarSummary(xAnom)
    printMinMax(xAnom, True)

    xAnom_sm   = calcDayAnomTLL (x, yyyyddd, xClmDay_sm)     
 xAnom_sm!0 = "time"
 xAnom_sm!2 = "lon"
 xAnom_sm!1 = "lat"

    xAnom_sm@long_name = "Anomalies from Smooth Daily Climatology"
    printVarSummary(xAnom_sm)
    printMinMax(xAnom_sm, True)

    delete( x )    ; no longer needed
dimx   = dimsizes(xAnom)
    ntim   = dimx(0)
    nlat   = dimx(1)
    mlon   = dimx(2)

;************************************************
; plotting parameters
;*******************************************G****
   if (PLOT) then
        LAT   = (/0/)
        LON   = (/0/)
        nPts  = dimsizes( LAT )

        plot  = new ( nPts, graphic)
        data2  = new ( (/2,366/), typeof(xClmDay), getFillValue(xClmDay))

        wks   = gsn_open_wks (wksType,wksName)    

        res                   = True                      ; plot mods desired
        res@gsnDraw           = False
        res@gsnFrame          = False
        res@trXMinF           =   1
        res@trXMaxF           = 366
        res@trYMinF           = 0
       ;res@tiMainString      = ""
       
        res@xyLineThicknesses = (/1.0, 2.0/)              ; make 2nd lines thicker
        res@xyLineColors      = (/"blue","red"/)          ; change line color
        res@xyMonoDashPattern = True                      ; all solid

        do np=0,nPts-1 
           data2(0,:) = xClmDay(:,0,0);{LAT(np)},{LON(np)})
           data2(1,:) = xClmDay_sm(:,0,0);{LAT(np)},{LON(np)})
           ;res@gsnCenterString = "lat="+LAT(np)+"  lon="+LON(np)
           res@tiXAxisString         = "Dia"
           res@tiYAxisString         = "%"
           plot(np)  = gsn_csm_y (wks,data2,res) ; create plot
        end do
       
        resP                     = True                ; modify the panel plot
        resP@txString            = vName+"-Climatologia "+yrStrt+"-"+yrLast
        resP@gsnMaximize         = True
        resP@gsnPaperOrientation = "portrait"
        gsn_panel(wks,plot,(/1,2/),resP)               ; now draw as one plot

       ;==========
       ; Plot anomalies for an arbitrarily selected near equatorial location
       ; Time: Oct 1, 1996 to April 1,1997  [arbitrary selection]
       ;==========
        LATX    = 0
        LONX    = 0

        yyyymmdd  = yyyyddd_to_yyyymmdd(yyyyddd)
      ;;yrfrac    = yyyymmdd_to_yyyyfrac (yyyymmdd, 0)
      ;;delete(yrfrac@long_name)

        xAnom@long_name    = "Anomalias Data Cruda (%)"   ; short labels for plot
        xAnom_sm@long_name = "Anomalias Filtro FFT (%)"

        ntBegin   = ind(yyyymmdd.eq.20050125)
        ntEnd     = ind(yyyymmdd.eq.20091005)
        
         monthLabels    = (/1,4,7,10/)
;        monNam = (/"Jan","Feb","Mar","Apr","May","Jun" \
;                  ,"Jul","Aug","Sep","Oct","Nov","Dec" /)
 monNam = (/"Ene","Feb","Mar","Abr","May","Jun" \
                  ,"Jul","Ago","Sep","Oct","Nov","Dic" /)

        xVal   = new(ntim, "integer", "No_FillValue") ; bigger than
        xLab   = new(ntim, "string", "No_FillValue")        ; needed
        xValm  = new(ntim, "integer" , "No_FillValue") ; bigger than

        ntm            = -1
        cr             = inttochar(10)                     ; carriage return
        do nt=ntBegin,ntEnd
         if (day(nt).eq.1) then
             ntm       = ntm + 1
             xVal(ntm) = xAnom&time(nt)
             xLab(ntm) = monNam(month(nt)-1)
             if (month(nt).eq.1) then
                 xLab(ntm) = xLab(ntm) + cr +sprinti("%0.4i", year(nt))
             end if
         end if
        end do

        rxy  = True      
        rxy@gsnDraw     = False
        rxy@gsnFrame    = False
        rxy@gsnYRefLine = 0.0                    ; create a reference line   
        rxy@gsnAboveYRefLineColor = "red"        ; above ref line fill red
        rxy@gsnBelowYRefLineColor = "blue"       ; below ref line fill blue
        rxy@xyLineThicknessF  = 2.0                                               
        rxy@vpHeightF  = 0.4                     ; resize            
        rxy@vpWidthF   = 0.8                  
        rxy@trXMinF           =   2005.0
        rxy@tmYRAutoPrecision=False
        rxy@tmXBPrecision=4
        ;rxy@tmXBMode   = "Explicit"
        ;rxy@tmXBValues = xVal(0:ntm)
        ;rxy@tmXBLabels = xLab(0:ntm)

      ;dia=ispan(0,dimsizes(ndvi)-1,1)
      ;ti=2005+(25+dia)/30
        fech=yyyymmdd_to_yyyyfrac(yyyymmdd,0.)
        plot2  = new ( 2, graphic)

        plot2(0)  = gsn_csm_xy (wks,fech\
                              ,xAnom(:,0,0),rxy) 
        plot2(1)  = gsn_csm_xy (wks,fech\
                              ,xAnom_sm(:,0,0),rxy) 
        resP@tiXAxisString         = "Dias desde 2005-01-25"
        resP@txString            = vName+"-Anomalias "
        resP@gsnMaximize         = True
        resP@gsnPaperOrientation = "portrait"

        gsn_panel(wks,plot2,(/2,1/),resP)   

    draw(plot2)
    end if
end
