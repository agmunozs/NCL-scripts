load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl" 
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl" 
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"

begin
  fname  = "PorcentajeLemnaPromedio.txt"
  data   = asciiread(fname,-1,"string")
  fecha  = stringtoint(str_get_cols(data,4,5))
  anomes = stringtoint(str_get_cols(data,0,5))
  ;ano    = stringtoint(str_get_cols(data,0,3))
  x0     = stringtofloat(str_get_cols(data,8,13))
  x0@_FillValue = -0.99

  Ndim   = dimsizes(x0)   ; number of elements 
  N      = Ndim(0)
  time1  = yyyymm_to_yyyyfrac(anomes, 0.5)
  time0  = yyyymm_time(2005, anomes(N-1), "integer")
  time   = time0({200501:anomes(N-1)})
  Ndim2  = dimsizes(time)   ; number of elements 
  N2     = Ndim2(0)

  fname2 = "ClorofilaA_Gio_NASA.txt"
  data2  = asciiread(fname2,-1,"string")
  clor   = new((/N2/),"float")
  cloro  = stringtofloat(str_get_cols(data2,10,15))
  clor(:N2-2)=cloro(:)
  clor(N2-1)=dim_avg(clor(:N2-2))
printVarSummary(clor)

  time2  = fspan(time1(0),time1(N-1),N2)

;print("Fechas inicial y final: "+fecha(0)+" - "+fecha(N-1))
print("Tiempos inicial y final:"+time(0)+ " - "+time(N2-1))

  prom = new((/N2/),"float")
  ic   = fecha(0)
  suma = 0. 
  k    = 1
  me   = 0

  do i=0,N-1
    if(fecha(i).eq.ic) then
       suma     = suma+x0(i) 
       k        = k+1
     else
       ic       = fecha(i)
       prom(me) = suma/k
       me       = me+1
       k        = 1
       suma     = x0(i)
    end if
    if(i.eq.(N-1)) then
       prom(me) = suma/k
    end if
  end do  

  prom@_FillValue = x0@_FillValue
  clor@_FillValue = x0@_FillValue

asciiwrite ("NDVI_medmens.dat", prom)

; Con el dtrend eliminamos tendencia y media:
  x  = prom ;dtrend_n(prom,False,0)
  xcl= clor ;dtrend_n(clor,False,0)

;Escritura en archivos ascii
  asciiwrite ("NDVI_ano_medmens.dat", time +"   "+x)
  asciiwrite ("clor_ano_medmens.dat", time +"   "+xcl)

;Calculo de correlacion y significancia estadistica
  spc = esccr( x, xcl,12)  ;spcorr
  siglvl= 0.05            ; a-priori specified sig level
  prob  = rtest(spc, N2, 0)
  print("Signif: "+(1.-prob))


;Entrenamiento del modelo mediante los parametros alfa y gamma:
  
; Requiramos derivada nula al inicio de la serie de tiempo.
;
ftsetp("sf1",1)
ftsetp("sl1", 0.0)
; Interpolacion con splines. Es una prueba para ver que tan bien funciona.
;
;yo = ftcurv(time2, prom, time2)

; Derivadas
;
lemnad = ftcurvd(time2, prom/prom(0), time2)
clorod = ftcurvd(time2, clor/clor(0), time2)

; Entrenamiento:
;parametros de eficiencia de captura y eficiencia de conversion reescalados:
;betab=1.
;deltab=1.
;alfa = new((/N2/),"float")
;gamma= new((/N2/),"float")
;alfa = lemnad/(prom/prom(0)) + betab*clor/clor(0)
;gamma= deltab*prom/prom(0)-clorod/(clor/clor(0))

; Entrenamiento:
;parametros de eficiencia de captura y eficiencia de conversion sin reescalar:
betab=1.
deltab=1.
;parametro de competencia (factor C1*N**2 en Lotka Volterra para competencia interespecifica):
C1=1e-2;
alfa = new((/N2/),"float")
gamma= new((/N2/),"float")
alfa = lemnad/(prom) + betab*clor+C1*prom
gamma= deltab*prom-clorod/(clor)

asciiwrite ("alfa.dat", alfa)
asciiwrite ("gamma.dat", gamma)

;Estacionales
alfa!0="time"
gamma!0="time"
estacl_2=month_to_season12(alfa(:71))
estacc_2=month_to_season12(gamma(:71)) 

printVarSummary(estacl_2)
;Calculo de correlacion y significancia estadistica
  spc2 = esccr( estacl_2, estacc_2,12)  ;spcorr
  siglvl= 0.05            ; a-priori specified sig level
  prob2  = rtest(spc2, N2, 0)
  print(1.-prob+"  "+spc2)
;********************************************************************************
; PLOTS
;********************************************************************************
 

  wks = gsn_open_wks("png","competencia")               ; open ps file
  wks2 = gsn_open_wks("png","corrlag") 
  wks3 = gsn_open_wks("png","deriv")

;**************Serie
  resI   = True
  resI@gsnCenterString     = "Competence. (Corr="+spc(0)+")"
  resI@gsnYRefLine   = 0
  ;res3@gsnAboveYRefLineColor = "Green"
  ;res3@gsnBelowYRefLineColor = "Blue"
  resI@xyLineColors   = (/"red","blue","green","yellow"/)
  resI@xyLineThicknesses = (/5.0,3.0,3.0,3.0/)
  resI@vpHeightF           = .4
  resI@vpWidthF            = .87
  resI@tiYAxisString = "NDVI Anomaly (% of surface)"
  resI@tiXAxisString = ""
  resI@tmXBTickStartF      = 2005
  resI@tmXBTickEndF        = 2012
  resI@tmXBPrecision       = 4
  resI@gsnMaximize = True
  resI@gsnPaperOrientation = "portrait"
  resI@trYMinF = 4
  resI@trYMaxF = -4

  resD = True
  resD@xyLineColors   = (/"blue","green","yellow"/)
  resD@xyLineThicknesses = (/5.0,3.0,3.0,3.0/)
  resD@tiYAxisString = "Clorophyll A Anomaly (mg m~S~-3~N~)"
  resD@gsnYRefLine   = 0
  resD@trYMinF = 20
  resD@trYMaxF = -20
 print(spc)

 
 res = True
  res@gsnCenterString     = "Correlation"
  res@tiXAxisString = "Lag (months)"
  plot = gsn_csm_xy2(wks,time2,x,xcl,resI,resD)
  plot2 = gsn_csm_y(wks2,spc,res)


  resI@gsnCenterString     = "Alfa (rojo) y Gamma (azul) --- Corr="+spc2(0)
  resI@gsnYRefLine   = 0
  resI@xyLineColors   = (/"red","blue","blue","blue"/);(/"red","red","blue","blue"/)
  resI@xyLineThicknesses = (/5.0,5.0,5.0,5.0/)
  resI@vpHeightF           = .4
  resI@vpWidthF            = .87
  resI@tiYAxisString = ""
  resI@tiXAxisString = ""
  resI@tmXBTickStartF      = 2005
  resI@tmXBTickEndF        = 2011
  resI@tmXBPrecision       = 4
  resI@gsnMaximize = True
  resI@trYMinF = 0  ;15
  resI@trYMaxF = 32 ;-10

printVarSummary(alfa)
printVarSummary(rmAnnCycle1D(alfa(:71)))

  plot3 = gsn_csm_xy(wks3,time2(:71),(/estacl_2,estacc_2/),resI)

;,prom/prom(0)-log(prom/prom(0))+alfa/gamma*(clor/clor(0)-log(clor/clor(0)))/),resI)
  

end
