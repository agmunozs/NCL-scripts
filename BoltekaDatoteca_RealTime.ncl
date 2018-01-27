;**************************************************************************************
;Centro de Modelado Cient’fico (CMC). La Universidad del Zulia   -- cmc.org.ve
;Observatorio Latinoamericano de Eventos Extraordinarios (OLE2)  -- ole2.org
;Copyleft: Angel G. Munoz S. -- agmunoz@cmc.org.ve
;Descripcion: script NCL del CMC para generar archivo netcdf con datos en reticula Boltek
;Journal: -This is Real Time version -- for RETRO see "BoltekaDatoteca.ncl"
;Version: 1.1 - May 4th, 2016. AGMS fixed: timestep problem, writing one NC per timestep
;!!!Esta versi—n borra todos los NETCDF presentes en el directorio cada vez que se ejecuta
;***************************************************************************************
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
begin

;************************************************
; Secci—n modificable por el usuario
;************************************************
  archi  		= "actual.dat"
  minlat		=8. ; lo obvio
  maxlat		=13.
  minlon		=-74.
  maxlon		=-70.
  reso			= 4.  ; Resolucion deseada en km
  detector		= "CMCMbo"

  ;###Archivos NetCDF salida####
  faltante 		= -999  ;bandera de valor faltante en archivos de salida
  diro 			= "./"   
  titulo	    = "Lightning Frequency (CG/IC)"
  sensor        = "Boltek_CMC_Mbo"
  autores       = "SIVIGILA (http://cmc.org.ve/sivigila.php -- Centro de Modelado Cientifico (CMC)) "

;------------------------------------------------
; FIN SECCIîN MODIFICABLE POR USUARIO
;------------------------------------------------

archsal	    = "_Boltek_"+detector+".nc"
archtemp="temp.dat"
system("rm -Rf *.nc")
system("rm -Rf "+archtemp)
system("cat "+archi+" | sed -e 's/,/ /g' > "+archtemp)

;************************************************
; Lectura archivo ASCII Boltek
;************************************************

  ncol   = numAsciiCol(archtemp)+4
  nrow   = numAsciiRow(archtemp)
  
  data = asciiread(archtemp,(/nrow,ncol/),"float")
  
  ano    = floattoint(data(0,0))
  mes    = floattoint(data(0,1))
  dia    = floattoint(data(0,2))
  hor    = floattoint(data(:,3))
  mit    = floattoint(data(:,4))
  seg    = floattoint(data(:,5))  
  lat    = data(:,6)
  lon    = data(:,7)
  ;pol    = data(:,8)				;no necesario por ahora
  tip   = floattoint(data(:,9))     ;0: CG --- 1: IC
  
;*****************************************************************************
; Definimos la ret’cula computacional
;*****************************************************************************
  lonr = fspan(minlon, maxlon,floattoint(abs(minlon-maxlon)*111.111/reso)) 
  nlon = dimsizes( lonr)
  latr = fspan(minlat, maxlat,floattoint(abs(maxlat-minlat)*111.111/reso))  
  nlat = dimsizes( latr)
  lonr!0="X"
  latr!0="Y"
  lonr@units    = "degree_east"
  latr@units    = "degree_north"
  
;*****************************************************************************
; Tiempos
;*****************************************************************************

  ntim 		= dimsizes(mit) 
  inter		= ispan(0, 1439, 10)   ; 144 timesteps, every 10 min
  ninter	= dimsizes(inter) 
  print(ntim)
  print(ninter)
  inter!0	= "T"
  inter@units="minutes since "+ano+"-"+sprinti("%0.2i",mes)+"-"+sprinti("%0.2i",dia)+" 00:00:00"

  t_evento	= new(ntim,double,faltante)
  
  do it=0,ntim-1
  	t_evento(it)	= cd_inv_calendar( ano, mes, dia, hor(it), mit(it), seg(it),  inter@units, 0 )
  end do
  ;print( t_evento(0:) )
  
  
;*****************************************************************************
; Campos
;*****************************************************************************
  cgf	= new((/ninter,nlat,nlon/),float,faltante)
  icf	= new((/ninter,nlat,nlon/),float,faltante)
  cgf	= 0
  cgf!0	="T"
  cgf!1	="Y"
  cgf!2	="X"
  cgf@units="frequency"
  icf	=cgf

  
  do it=0,ntim-1
  	indti       = closest_val(t_evento(it),inter)  				;valor m‡s cercano en tiempo
  	if (inter(indti).ne.0.and.t_evento(it)/inter(indti).lt.1.) then     ;revisamos si pertenece a 10 min anteriores
  		indti 	= indti-1
  	end if
  	indla       = closest_val(lat(it),latr)						;valor m‡s cercano en lat
  	indlo       = closest_val(lon(it),lonr)						;valor m‡s cercano en lon
  	if (tip(it).eq.0) then
  	 ;if (cgf(indti,indla,indlo).eq.faltante) then        
  	 ;  cgf(indti,indla,indlo)=0							;if detection happened, it isn't a missing
  	 ;end if
  	  cgf(indti,indla,indlo)=cgf(indti,indla,indlo)+1.
  	else
  	 ;if (icf(indti,indla,indlo).eq.faltante) then
  	 ;  icf(indti,indla,indlo)=0							;if detection happened, it isn't a missing
  	 ;end if
  	  icf(indti,indla,indlo)=icf(indti,indla,indlo)+1.
  	end if
  end do
  
  last=indti
;*****************************************************************************
; Escritura archivo NETCDF
;*****************************************************************************

   do it=0,last
   	it1=it+1                    
	system("/bin/rm -f " + diro + ano+""+sprinti("%0.2i",mes)+""+sprinti("%0.2i",dia)+"_"+sprinti("%0.3i",it1)+archsal) ; borrar si existe (s—lo ese d’a)   
	fout  = addfile (diro +ano+""+sprinti("%0.2i",mes)+""+sprinti("%0.2i",dia)+"_"+sprinti("%0.3i",it1)+archsal, "c")  ; abrir archivo de salida

    setfileoption(fout,"DefineMode",True)

    fAtt               = True            
	
	fAtt@source		   = sensor
	fAtt@authors	   = autores 
	fAtt@title         = titulo  
	fAtt@creation_date = systemfunc ("date")        
	fileattdef( fout, fAtt )            ; copy file attributes    

    dimNames = (/"T", "Y", "X"/)  
	dimSizes = (/ -1   ,  nlat,  nlon/) 
	dimUnlim = (/ True , False, False/)   
	filedimdef(fout,dimNames,dimSizes,dimUnlim)

       filevardef(fout, "T"  ,"integer",(/"T"/) ) 
       filevardef(fout, "Y"  ,typeof(latr),getvardims(latr))                          
       filevardef(fout, "X"  ,typeof(lonr),getvardims(lonr))                          
       filevardef(fout, "CG" ,typeof(cgf),(/  "Y", "X" /) ) ;getvardims(cgf))  
       filevardef(fout, "IC" ,typeof(icf),(/  "Y", "X" /) ) ;getvardims(icf))
    ;===================================================================
    ; Copy attributes associated with each variable to the file
    ; All attributes associated with each variable will be copied.
    ;====================================================================
       filevarattdef(fout,"T" ,inter)                     
       filevarattdef(fout,"Y"  ,latr)                     
       filevarattdef(fout,"X"  ,lonr)                     
       ;filevarattdef(fout,"CG",cgf)                      
       ;filevarattdef(fout,"IC",icf)
       

       fout->T   	= (/inter(it)/)     
       fout->Y   	= (/latr/)
       fout->X    	= (/lonr/) 
       fout->CG  	= (/cgf(it,:,:)/)
       fout->IC  	= (/icf(it,:,:)/)

   end do
  
 end

