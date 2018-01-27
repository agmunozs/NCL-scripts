;****************************************************
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRFUserARW.ncl"
;****************************************************
begin
;***************************
; read in data
;***************************

a = addfile("1998_2008.nc","r")
b = addfile("CATA_TempPrec.nc","r")
c = addfile("CATA_AdHum.nc","r")
d = addfile("CATA_radiacion.nc","r")
f= addfile ("Cata1998_CAPEvars.nc","r")
g = addfile("CATA_pres.nc","r")
h = addfile("LISOTD_LRTS_V2.2.nc","r")

 anoa=1998
 mesa=1
 lonmin=-73.
 latmin=8.
 lonmax=-70
 latmax=11.
 llres = True
  llres@ReturnInt = False   ; Return integer values
  locmin = wrf_user_ll_to_ij(f, lonmin, latmin, llres)
  locmax = wrf_user_ll_to_ij(f, lonmax, latmax, llres)
  locXmin = locmin(0)
  locYmin = locmin(1)
  locXmax = locmax(0)
  locYmax = locmax(1)

 XLAT  = a->XLAT
 XLON  = a->XLONG
 EPDC  = a->CAPECIN
    U  = c->Umen
    V  = c->Vmen
  Tem  = b->T
 PSFC  = g->PSFCmen
  PSFC = PSFC / 100.
  BAL  = d->DifRadmen
  flash = h->LRTS_COM_FR ;LRTS_LIS_FR
 x     = h->Day_since_Jan_1_95_COM  ;Day_since_Jan_1_95_LIS
 lat   = h->Latitude
 lon   = h->Longitude  
           
time = fspan(1998,2009,132)
timer= fspan(1998,2006,114)
 tini=1096;  Para 98=1096; Aqui va el numero inicial de dias
 anoi=1998; aqui t final
 superf=277.*277.;superficie de la celda en km2
  suma=(greg2jul(1998,1,1,-1))
x2=anoi+(x-tini)/365

medeCAPE = new(dimsizes(EPDC(:,0,0,0,0)),"float")
do i=0,dimsizes(EPDC(:,0,0,0,0))-1
   medeCAPE(i)=dim_avg_Wrap(dim_avg_Wrap(EPDC(i,0,0,locYmin:locYmax,locXmin:locXmax)))
end do


anoa=1998
mesa=1

 ndim = dimsizes(flash)
print(ndim)
  ntim=ndim(2)
  nlat=ndim(0)
  nlon=ndim(1)
  nmos = 12.
  ntimes=stringtoint((ntim-1096)/365.*nmos)
print(ntimes)
FLASHmen = new((/nlat,nlon,ntimes/),float)

  ini=0
  anoi=anoa
  mesi=mesa
  print("Procesando Valores Medios")
  do i=0,ntimes-1
    dt=days_in_month(anoi,mesi)
    FLASHmen(:,:,i)=dim_avg_Wrap(flash(Latitude|:,Longitude|:,Day_since_Jan_1_95_COM|ini:ini+dt))
  ini=ini+dt
print("Mes:"+mesi+"   Año:"+anoi)
    mesi=mesi+1
    if(mesi.eq.13) then
      mesi=1
      anoi=anoi+1
      print("Mes:"+mesi+"   Año:"+anoi)
    end if
  end do

medeFLASH = new(dimsizes(flash(0,0,:)),"float")
do i=1096,dimsizes(flash(0,0,1096:))-1
   medeFLASH(i)=dim_avg_Wrap(dim_avg_Wrap(FLASHmen(0,0,:)))
end do



;***************************
; plot parameters
;***************************
  wks = gsn_open_wks("ps","SerieFlash")
  wks1 = gsn_open_wks("ps","SerieFlash1")
; left variable
  resL           = True

  resL@gsnFrame  = False           ; Don't draw plot or advance frame
  resL@gsnDraw   = False

  resL@vpXF      = 0.15
  resL@vpWidthF  = 0.58            ; Set width and height of plot.
  resL@vpHeightF = 0.2

  resL@tmXBMode            = "Manual"
  resL@tmXBTickStartF = 1998.
  resL@tmXBTickEndF = 2009.
  resL@trXMinF = resL@tmXBTickStartF
  resL@trXMaxF = resL@tmXBTickEndF

  ;resL@trXMinF = 1998.                 ; Set min/max of X axis.
  ;resL@trXMaxF = 2009.

  resL@tmXBValues        = time

  resL@xyLineColor       = "NavyBlue"   ; To distinguish the lines
  resL@xyLineThicknesses = 2.      ; thicker line

; tickmarks:
  resL@tmXTOn              = False   ; Turn off top tickmarks

  resL@tmYLMinorOn         = False   ; Turn off left minor tickmarks

  resL@tmXBLabelFontHeightF = 0.01

  resL@tmXMajorGrid                = True
  resL@tmYMajorGrid                = True
  resL@tmXMajorGridThicknessF      = 1.0
  resL@tmYMajorGridThicknessF      = 1.0
  resL@tmXMajorGridLineDashPattern = 2

 resL@tiMainString                = "Time Series of Mean CAPE and RB"
  resL@tiYAxisString               = "CAPE (J/kg)"  ; axis string
  resL@tiYAxisFontColor            = resL@xyLineColor
  resL@tiYAxisFontHeightF          = 0.01

; right variable
  resR                        = True
  resR@xyLineColor            = "Brown"
  resR@xyDashPatterns         = 1                   ; dashed line for 2nd
  resR@xyLineThicknesses      = 2                   ; thicker line

  resR@trXMinF = resL@trXMinF
  resR@trXMaxF = resL@trXMaxF

  resR@tiYAxisString    = "Zonal Winds (m/s)"   ; axis string
  resR@tiYAxisFontColor = resR@xyLineColor
  resR@trYReverse              = True
  resR@tmYRMinorOn         = False   ; Turn off right minor tickmarks

; Create the plot with the three variables and three different axes.
; Nothing is drawn at this point.
;
  plot= gsn_csm_xy2(wks,time,medeCAPE,medeFLASH,resL,resR)
  plot1= gsn_csm_xy(wks1,timer,medeFLASH,resR)

draw(plot)
draw(plot1)
end



