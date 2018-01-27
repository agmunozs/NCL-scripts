; *********************************************************
; NCL Script for reading a GrADS file and writing 
; in CPTv>=10 format
; Angel G. Munoz (agmunoz@iri.columbia.edu)
; International Research Institute for Climate and
; Society (IRI). Columbia University
; AND
; Observatorio Latinoamericano de Eventos Extraordinarios
; (OLE2) 
; Dec 2011
; *********************************************************
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl" 
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

begin
;**********************************************************
; Begin of user modification section
;**********************************************************
  diri  = "./"
  fili  = "psa.jjas"                      ; DSET
  fName = diri+fili       
  filen= "psa.jjas"
  vari="t"
  varname = "precip"         
  units= "Celsius"                 
  timty= 2   ; 1: monthly (all), 2: only one month/year, 3: seasonal (not implemented)
  tispe= 7 ; speficy month as integer (e.g. 1 for Jan, if timty = 2)
  yrini= 1979
  yrend= 2010
  missing=-999000000.00
  
  PLOT  = False     ; True if plots are desired
  NC    = False		; True if netcdf is desired
  
;INFO from the .ctl file:
  nlat  = 72                            ; YDEF
  mlon  = 144                              ; XDEF
  lonmin = 1.25
  latmin = -88.75
  xstep = 2.5
  ystep = 2.5
  
  year  = yrini                            ; TDEF
  ntim  = 32                               ; time steps  yrs
  ;nmos  = 12
  
;**********************************************************
; End of user modification section
;**********************************************************
                                          ; not required
  ;time  = new (ntim*365, float  ,"No_FillValue") ; generate a "time" variable
  ;date  = new (ntim*365, integer,"No_FillValue") ; generate a "date" variable

  time0=yyyymm_time(yrini,yrend,"integer")  
  if(timty.eq.1) then
  time=time0
  delete(time0)
  end if
  if(timty.eq.2) then
  time=time0(tispe-1::12)
  end if
  print(time)
  yearo= str_get_cols(time, 0, 3)
  mont= str_get_cols(time, 4, 5)


  time!0         = "time"
  time@long_name = "time"
  time@units     = "???"           ; "yyyy.fraction_of_year"
  time&time      = time
                                          ; generate lat/lon
  lon = ispan(0,mlon-1,1)*xstep + lonmin
  lon!0 = "lon"
  lon@long_name = "longitude"
  lon@units     = "degrees_east"

  lat = ispan(0,nlat-1,1)*ystep + latmin 
  lat!0 = "lat"
  lat@long_name = "latitude"
  lat@units     = "degrees_north"
                                          ; create an array to contain data
  UNDEF = missing                           ; UNDEF
  x     = new ( (/ntim,nlat,mlon/), float, UNDEF)
  x!0   = "time"
  x!1   = "lat"
  x!2   = "lon"
  x&time=  time
  x&lat =  lat
  x&lon =  lon

  x@long_name = varname ; VARS
  x@units     = units
                                    
 ;setfileoption("bin","ReadByteOrder","Native")       ; a033 default
 setfileoption("bin","ReadByteOrder","LittleEndian")
 ; setfileoption("bin","ReadByteOrder","BigEndian")
                                         
                                          
;setfileoption("bin","RecordMarkerSize",8)
;setfileoption("bin","ReadByteOrder","BigEndian")
;print(fbinnumrec(fName))

  do nt=0,ntim-1                         ; the ::-1 changes the latitude order
     x(nt,:,:) = fbindirread(fName, nt, (/nlat,mlon/), "float")
     ;print(fbindirread(fName, nt, (/nlat,mlon/), "float"))
  end do

  printVarSummary(x)
  print ("min(x)="+min(x))
  print ("max(x)="+max(x))

  if (NC) then
      nline= inttochar(10)

      diro = "./"
      filo = filen + ".nc"
      system ("/bin/rm -f "+diro+filo)  ; remove any pre-existing file

      ncdf = addfile (diro+filo, "c")
     ;setfileoption(ncdf,"DefineMode",True)               ; a033 [most efficient]
      
      globeAtt         = 1      ; global [file] attributes
      globeAtt@title   = "Precipitation"
      globeAtt@source  = fili 

      globeAtt@story   = "GrADS to CPT -- AGMunoz (agmunoz@iri.columbia.edu)"
  
      globeAtt@creation_date= systemfunc ("date" ) 

      fileattdef( ncdf, globeAtt )
                                    
      dimNames = (/"time", "lat", "lon" /)  
      dimSizes = (/  -1  ,  nlat,  mlon /) 
      dimUnlim = (/ True , False, False /)   
      filedimdef(ncdf, dimNames  , dimSizes,  dimUnlim )

      filevardef   (ncdf, "time"  , typeof(time) , "time" )
      filevarattdef(ncdf, "time", time)
     
      filevardef   (ncdf, "lat", "float", "lat")
      filevarattdef(ncdf, "lat", lat)
     
      filevardef   (ncdf, "lon", "float", "lon")
      filevarattdef(ncdf, "lon", lon)
      
      filevardef(ncdf, "var"  , typeof(x) , (/"time", "lat" , "lon"/) )
      filevarattdef(ncdf, "var", x)

      ncdf->time   = (/ time /)
      ncdf->lat    = (/ lat /)
      ncdf->lon    = (/ lon /)
      ncdf->var   = (/ x /) 
  end if

  if (PLOT) then
      wks    = gsn_open_wks("ps","GPCP_grads")      ; open a ps file
      colors = (/"white","black"             \    
                ,"azure1","beige","lavender" \
                ,"PaleGreen","SeaGreen3","LightYellow" \
                ,"Yellow","HotPink","Red", "Purple"/)      
      gsn_define_colormap(wks, colors)              ; generate new color map
    
      res = True
      res@gsnMaximize          = True               ; make large
      res@cnFillOn             = True               ; turn on color fill
      res@cnLinesOn            = False              ; turn off contour lines
      res@cnLineLabelsOn       = False              ; turn off line labels
      res@cnLevelSelectionMode = "ExplicitLevels"   ; set explicit contour levels
      res@cnLevels             = (/ 0.5, 1.0, 3.0 \ ; set unequal contour levels
                                  , 5.0, 7.0,10.0 \
                                  ,15.0,25.0,50.0 /)
      res@mpFillOn             = False              ; turn off gray continents
      res@mpCenterLonF         = 180                ; Centers the plot at 180 
    
      x = mask(x, x.eq.0. , False)            
     
      do nt=0,1    ; ntim-1
         res@gsnCenterString      = date(nt)
         plot = gsn_csm_contour_map_ce(wks,x(nt,:,:),res)  
      end do

  end if

;**********************************************************

wlon=new((/1/), "string")
wlat=new(nlat, "string")
v   = x
v__FillValue = missing
delete(x)

;**********************************************************
; Write data
;**********************************************************

system("rm -Rf "+filen+".txt")

  asciiwrite("algo", "xmlns:cpt=http://iri.columbia.edu/CPT/v10/")
  system("cat algo >>"+filen+".txt")
  asciiwrite("algo", "cpt:nfields=1")
  system("cat algo >>"+filen+".txt")
  
  ltim="cpt:T		"
  do it=0,ntim-1
  	ltim=ltim+yearo(it)+"   "
  end do
asciiwrite("algo",ltim)
  system("cat algo >>"+filen+".txt")

do it=0,ntim-1
  print("processing: "+yearo(it)+ " " +mont(it))

	asciiwrite("algo", "cpt:field="+vari+", cpt:T="+yearo(it)+"-"+mont(it)+", cpt:nrow="+nlat+", cpt:ncol="+mlon+", cpt:row=Y, cpt:col=X, cpt:units="+units+", cpt:missing="+v@_FillValue)
	system("cat algo >>"+filen+".txt")
	wlon = "        "  
	do nl=0,mlon-1        
            wlon = wlon + sprintf("%7.2f   ",lon(nl))
 	end do
 	asciiwrite("algo",wlon)
	system("cat algo >>"+filen+".txt")


  do nl=0,nlat-1        
       wlat(nl) = sprintf("%7.1f ",lat(nl))
    do ml=0,mlon-1
	 wlat(nl) = wlat(nl) + sprintf("%7.2f   ", v(it,nl,ml))
    end do
  end do
 	asciiwrite("algo",wlat)
	system("cat algo >>"+filen+".txt")

end do

system("rm -Rf algo")
end
