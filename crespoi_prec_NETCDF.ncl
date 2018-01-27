;**************************************************************************************
;Centro de Modelado Científico (CMC). La Universidad del Zulia   -- cmc.org.ve
;Observatorio Latinoamericano de Eventos Extraordinarios (OLE2)  -- ole2.org
;Copyleft: Angel G. Munoz S. -- agmunoz@cmc.org.ve
;Descripcion: script NCL del OLE2 y CMC para generar archivo netcdf con datos en reticula
;***************************************************************************************
;---------------------------------------------------------------------------------------
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"
load "./shapefile_mask_data.ncl"
;--------------------------------------------------------------------------------------- 

begin
  archivo = "data_PAN_P.dat"  ;Nombre del archivo de datos (entrada)
  npts = 84 ; Num estaciones
  ntim = 408; Num meses a procesar
  faltante = -999.99  ;bandera de valor faltante en archivo de entrada
  minlat=5 ; lo obvio
  maxlat=11
  minlon=-86
  maxlon=-75
  anoini=1981
  reso= 25 ; Resolucion deseada en km
  relleno = 0; 0: sin relleno de Poisson; 1: relleno con Poisson
  rscan = (/0.5,0.3,0.1/)   ;Radios de influencia
  exper = 1 ; 1 para indicar que el archivo NetCDF de salida 
  			; incluya informaci—n del experimento en el nombre
  ;###Archivo NetCDF####
  archsal       = "climatologia_mens_p.nc"
  titulo        = "Precipitacion Mensual (mm)"
  metodo        = "Cressman (1959) con Poisson (relajacion). Resolucion "+reso+" km"
  instituto     = "ALGO -- Centro de Modelado Cientifico (CMC) - LUZ"
  autores       = "Angel G. Munoz (agmunoz@cmc.org.ve) "
  creacion      = systemfunc ("date")


;*****************************************************************************
;Seccion Lectura Archivo ASCII
;*****************************************************************************
  data = asciiread(archivo,(/ntim*npts,6/),"float")

; Años y meses
  ano1d= data(::npts,0)  ;con el ::ALGO hacemos que lea cada ALGO pasos
  mes1d= data(::npts,1)

; Para las latidudes de las estaciones, basta adquirir las primeras ntim lineas:
  lat1d = data(0:npts-1,3)

; Idem para las longitudes
  lon1d  = data(0:npts-1,4)

;***************************************************************************
;Abrimos archivo CRU para first guess del analisis objetivo
;**************************************************************************
  f    = addfile ("precip_CRU_mes_71-00.nc", "r")
  pr   = f->prcp(36:,:,:)
  TCRU = f->T(36:)
  lati = f->Y
  loni = f->X
  pr!0="time"
  pr!1="lat"
  pr!2="lon"
;printVarSummary(pr)

;*****************************************************************************
; Definimos la nueva malla para el Cressman
;*****************************************************************************
  lon = fspan(minlon, maxlon,floattoint(abs(minlon-maxlon)*111.111/reso)) ;Creamos la malla (longitud, res=30 km)
  nlon = dimsizes( lon)
  ;print(lon)
  lat = fspan(minlat, maxlat,floattoint(abs(maxlat-minlat)*111.111/reso))  ;Creamos la malla (latitud, res=30 km)
  nlat = dimsizes( lat)
 ; printVarSummary(lat)
 ; printVarSummary(lon)
  X=lon
  Y=lat
;****************************************************
;Interpolacion a la malla del Cressman de la data CRU
;****************************************************
pre=linint2 (pr&lon, pr&lat, pr, False, lon, lat, 0)
pre@_FillValue=faltante
pre@lon=lon
pre@lat=lat
pre!0="time"
pre!1="lat"
pre!2="lon"

;printVarSummary(pre)
;print(pre)

;*****************************************************************************
; Definimos unos arreglos a usar
 temp1D=new((/ntim,npts/),typeof(data))
 temp2D=new((/ntim,npts,npts/),typeof(temp1D))

; temp1D es arreglo de (tiempo,estacion)!!!!

; Ciclo de lectura y almacenamiento en arreglos para la variable física 
  do i=0,ntim-1 
   temp1D(i,:)   = data(i*npts:(i+1)*npts-1,5)            
;   temp2D(i,:,:) = onedtond(temp1D(i,:),(/npts,npts/)) ; convertimos a 2D
  end do
;print(temp1D)
;**************************************************
;Esto Define los valores faltantes (superfluo)
  temp1D@_FillValue    = faltante

;**************************************************
; Calculo de la anomalia estandarizada
;**************************************************
; Definimos arreglos a usar
  climo=new((/12,npts/),typeof(data))
  de   =new((/12,npts/),typeof(data))

;printVarSummary(temp1D)
;print(temp1D(:ntim-1:12,0))
  do i=0,npts-1
   do mes=0,11
     climo(mes,i) = dim_avg_n_Wrap(temp1D(mes::12,i),0)
     de   (mes,i) = dim_stddev_n_Wrap(temp1D(mes::12,i),0)
   end do
  end do
;printVarSummary(climo)

;**************************************************
; Seccion Cressman. Ajustar los radios de influencia al inicio (rscan)
;**************************************************
  grid =new((/nlat,nlon,ntim/),typeof(temp1D))
  T=ispan(0,ntim-1,1)
  T@units ="months since "+anoini+"-01-01"  ;inicio de mes para representar el mes 
; T@units ="months since "+anoini+"-01-15"  ;mitad de mes para representar el mes
  T@description="time"
  ;T@calendar = "365_day"
  grid!0="lat"
  grid!1="lon"
  grid!2="T"
  grid&T=T
  grid&lon=lon
  grid&lat=lat
  grid@_FillValue    = faltante
  printVarSummary(grid)

 ;Para sd y climo:
  km=0
  do nmo =0,ntim-1
    if(nmo-km.ge.12) then
      km=km+12
    end if
    zVal = temp1D(nmo,:)
    opt = True
    ;opt@guess = pre(nmo,:,:)  ; CRU como first guess, pasamos a mm
    opt@zonal = False ; esta opción (recomendada) permite que FG = medias zonales
    grid(:,:,nmo) = obj_anal_ic_deprecated(lon1d,lat1d,zVal,X,Y, rscan, opt) 
  end do
  grid_cl=clmMonLLT(grid)   
  grid_sd=stdMonLLT(grid)

 ;Calculamos la estandarizacion de la variable
  km=0
  do nmo =0,ntim-1
    if(nmo-km.ge.12) then
      km=km+12
    end if
    zVal = (temp1D(nmo,:)-climo(nmo-km,:))/de(nmo-km,:)
    opt = True
    ;opt@guess = pre(nmo,:,:)  ; CRU como first guess, pasamos a mm
    opt@zonal = False ; esta opción (recomendada) permite que FG = medias zonales
    grid(:,:,nmo) = obj_anal_ic_deprecated(lon1d,lat1d,zVal,X,Y, rscan, opt) 
  end do
	

 ;Calculamos la variable original en la reticula
  km=0
  do nmo =0,ntim-1
    if(nmo-km.ge.12) then
      km=km+12
    end if
    grid(:,:,nmo)=grid(:,:,nmo)*grid_sd(:,:,nmo-km)+grid_cl(:,:,nmo-km)
  end do

;No se aceptan precipitaciones negativas creadas por el proceso de interpolacion.
  grid=grid>0

  grid&lon@units    = "degree_east"
  grid&lat@units     = "degree_north"

;**************************************************
; Atributos requeridos para evitar warnings fastidiosos
;**************************************************
  ;pre&lon@units    = "degree_east"
  ;pre&lat@units     = "degree_north"

;******************************
;Ahora rellenamos con Poisson relajado
;********************************
if(relleno.eq.1) then
  print("Resolviendo Ec. de Poisson via relajacion")

  guess     = 1                ; use zonal means
  is_cyclic = False             ; cyclic [global]
  nscan     = 1500             ; usually much less than this
  eps       = 1.e-2            ; variable dependent
  relc      = 0.6              ; relaxation coefficient
  opt2       = 0                ; not used

 do nme=0,ntim-1
   poisson_grid_fill(grid(:,:,nme), is_cyclic, guess, nscan, eps, relc, opt2)
 end do
end if

;***********************************
;Mascaras
;***********************************
 ;b = addfile("landsea.nc","r")   ; read in land sea mask basemap file
 ;b = addfile("SampleLandSeaMask.nc","r")
 ;masc=b->LSMASK
 ;mascara  = landsea_mask(masc,lat,lon)
 ;printVarSummary(grid)
 ;printVarSummary(mascara)
 ; Create a mask array the same size as "u", using
; lat/lon data read off a shapefile.
;
  shpfile   = "GSHHS_h_L1.shp"
  opt             = True
  opt@return_mask = True

;  land_mask = shapefile_mask_data(grid(:,:,0),shpfile,opt)
  
; do ti=0,ntim-1
;    ;grid(:,:,ti) = mask(grid(:,:,ti),mascara.eq.0,False)
;    grid(:,:,ti)  = where(land_mask.eq.1,grid(:,:,ti),faltante)
; end do

 grid@_FillValue    = faltante
 grid@missing_value = grid@_FillValue
;**************************************************************************
if(exper.eq.1) then
        filo = rscan(0)+"-"+rscan(1)+"-"+rscan(2)+"_"+archsal             ; Archivo salida
        if(relleno.eq.1) then
        	filo = rscan(0)+"-"+rscan(1)+"-"+rscan(2)+"_poisson_"+archsal
        end if
   else
        filo = archsal      
end if
        system("/bin/rm -f " + filo)    ; remove if exists
        fout1  = addfile (filo, "c")  ; open output file

    ;===================================================================
    ; explicitly declare file definition mode. Improve efficiency.
    ;===================================================================
        setfileoption(fout1,"DefineMode",True)

    ;===================================================================
    ; create global attributes of the file
    ;===================================================================
        fAtt               = True            ; assign file attributes
		fAtt@creacion      = systemfunc ("date")
        fAtt@autores       = autores
        fAtt@instituto     = instituto
        fAtt@radiosinflue  = rscan(0)+"-"+rscan(1)+"-"+rscan(2)
        fAtt@metodo        = metodo
        fAtt@titulo        = titulo
        

       ;filedimdef(fout1,"time",-1,True)
       fileattdef( fout1, fAtt )            ; copy file attributes

       ;fout1->T        = T
       fout1->precip   = grid

;*********************************************************************************
  print("Fin exitoso del proceso - Ha sido un gusto para CMC servirle.")
;*********************************************************************************

end
