; *********************************************************
; NCL Script for reading a 1-station daily file and writing 
; in CPTv>=10 format
; Angel G. Munoz (agmunoz@iri.columbia.edu)
; International Research Institute for Climate and
; Society (IRI). Columbia University
; AND
; Observatorio Latinoamericano de Eventos Extraordinarios
; (OLE2) 
; May 2015
; *********************************************************
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

begin
;**********************************************************
; Begin of user modification section
;**********************************************************
filein="Barbados_CIMH_wet_days.txt"
fileou="CPT_freq_mon_ssn.txt"
var  = 1   ; 1: precip, 2 : temp
units= "counts"
missing = -999.
ystart=1969
yend  =2014
npts = 1 ; Num estaciones
ntim = 12418; Num total pasos de tiempo (num líneas del archivo estación en formato Climdex)
archlis=asciiread("lista.ls",(/npts,1/),"string")
;**********************************************************
; End of user modification section
;**********************************************************
  err = NhlGetErrorObjectId()
  setvalues err
  "errLevel" : "Fatal"    ;only report Fatal errors
  end setvalues
;**********************************************************
; Read in data
;**********************************************************

  dir="test"
  data0 = asciiread(dir + "/" + filein,(/ntim,4/),"float") ;4 columns !!!


  lat1d = data(0:npts-1,0)
  lon1d  = data(0:npts-1,1)

  col=3  ;for first parameter in file
  varl=new((/npts,ntim,4/),"float") ;4 columns !!!  change to 6 if original Climdex format
  do i=0,npts-1 
   algo = asciiread(dir + "/" + archlis(i,0),(/ntim,4/),"float") ;4 columns !!!
   varl(i,:,:) = algo
  end do

  var=varl(:,:,col)
  var!0="station"
  var!1="time"

  var@_FillValue    = missing

;Esto es para definir tiempos
 yyyymmdd  = yyyymmdd_time(ystart, yend, "integer")
  ntim  = dimsizes(yyyymmdd)
  yyyy  = yyyymmdd/10000 
  mmdd  = yyyymmdd-yyyy*10000                ; mmdd = yyyymmdd%10000
  mm    = mmdd/100
  dd    = mmdd-mm*100                        ; dd   = mmdd%100
  hh    = dd                                 ; create arrays [*] of required size
  mn    = dd
  sc    = dd
  hh    = 0                                  ; array syntax
  mn    = 0
  sc    = 0
  var&time = cd_inv_calendar(yyyy,mm,dd,hh,mn,sc,"days since 1961-01-01 00:00", 0)
  delete(yyyymmdd)
  delete(yyyy)
  delete(mmdd)
  delete(mm)
  delete(dd)
  delete(hh)
  delete(mn)
  delete(sc)
 printVarSummary(var)
;**************************************************
;Cálculo de valores mensuales
 yyyymm = yyyymm_time(ystart, yend, "integer")
 yyyy  = yyyymm/100
 mm  = yyyymm-yyyy*100               ; mmdd = yyyymmdd%10000
 dd    = mm
 dd    = 15
 hh    = dd                                 ; create arrays [*] of required size
 mn    = dd
 sc    = dd
 hh    = 0                                  ; array syntax
 mn    = 0
 sc    = 0
 
 nmon=(yend-ystart+1)*12
 vmo=new((/npts,nmon/),"float")
 do i=0,npts-1 
    vmo(i,:)= calculate_monthly_values(var(i,:), "sum", 0, False)
 end do
 vmo!0="station"
 vmo!1="time"
 vmo&time = cd_inv_calendar(yyyy,mm,dd,hh,mn,sc,"days since 1961-01-01 00:00", 0)
printVarSummary(vmo)
;**********************************************************
; Write out data
;**********************************************************

  system("rm -Rf "+filen+"_"+seasons(sea)+"_stn.txt")
  asciiwrite("algo", "xmlns:cpt=http://iri.columbia.edu/CPT/v10/")
  system("cat algo >>"+filen+"_"+seasons(sea)+"_stn.txt")
  asciiwrite("algo", "cpt:nfields=1")
  system("cat algo >>"+filen+"_"+seasons(sea)+"_stn.txt")

	asciiwrite("algo", "cpt:field="+vari+", cpt:nrow="+ntim+", cpt:ncol="+nstat+", cpt:row=T, cpt:col=station, cpt:units="+units+", cpt:missing="+v@_FillValue)
	system("cat algo >>"+filen+"_"+seasons(sea)+"_stn.txt")

	wlon = "        "  
	do nl=0,nstat-1        
            wlon = wlon + "     "+stID(nl)
 	end do
 	asciiwrite("algo",wlon)
	system("cat algo >>"+filen+"_"+seasons(sea)+"_stn.txt")

	wlon = "cpt:X  "  
	do nl=0,nstat-1        
            wlon = wlon + "     "+loni(nl)
 	end do
 	asciiwrite("algo",wlon)
	system("cat algo >>"+filen+"_"+seasons(sea)+"_stn.txt")

	wlon = "cpt:Y  "  
	do nl=0,nstat-1        
            wlon = wlon + "     "+lati(nl)
 	end do
 	asciiwrite("algo",wlon)
	system("cat algo >>"+filen+"_"+seasons(sea)+"_stn.txt")

    do i=sea,ntim-3,12
	moi=sea+1
	mof=sea+3
	if(mof.gt.12) then
	  mof=mof-12
	end if
	wlon = "    "+year(i)+"-"+sprinti("%2.2i",moi)+"/"+sprinti("%2.2i",mof)
	do nl=0,nstat-1      
	    ;falgo=dim_sum_n(vi(i:i+2,nl),0)
	;if second method (fill only the missing values with Cressman):
	     falgo=dim_sum_n(temp1D(i:i+2,nl),0)
	    if(falgo.lt.0) then 
		falgo=temp1D@_FillValue
	    end if    
               wlon = wlon + "     "+falgo
 	end do
 	asciiwrite("algo",wlon)
	system("cat algo >>"+filen+"_"+seasons(sea)+"_stn.txt")
    end do

system("rm -Rf algo")
end
