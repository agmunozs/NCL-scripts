load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
begin
;***********************************************
; get variable names from nc files
;***********************************************
   f1       = addfile("Data.D/sa23.daily.1.1940-2011.nc","r")   
   f2       = addfile("Data.D/tn_61-00.nc","r")
   f3       = addfile("Data.D/tx_61-00.nc","r")

;***********************************************
;Interpolacion
;***********************************************
; create output netcdf file
;*********************************************** 
   system("rm SESA_prctntx_1961_2000_1deg.nc") ; remove any pre-existing file
   ncdf_out = addfile("SESA_prctntx_1961_2000_1deg.nc" ,"c")       ; create output netCDF file


     precip=f1->precip(7671:22280,{-40.25:-19.75},{-70:-45})
     precip@_FillValue=-9999.
     precip@missing_value=precip@_FillValue
     precip@Origin="Brant Liebman precipitation dataset for South America"
     precip@scale_factor=1.0
     lat1  =f1->lat({-40.25:-19.75})
     lon1  =f1->lon({-70:-45})
     tn    =f2->tn
     tx    =f3->tx
     tn@missing_value=precip@_FillValue
     tx@missing_value=precip@_FillValue
     tn@Origin="Barbara Tencer minimum temperature dataset"
     tx@Origin="Barbara Tencer maximum temperature dataset"
     ncdf_out->time  = f2->time
     ncdf_out->prec = precip
     ncdf_out->tmi  = area_hi2lores_Wrap(tn&longitude,tn&latitude,tn,False,1,lon1,lat1,False)
     ncdf_out->tma  = area_hi2lores_Wrap(tx&longitude,tx&latitude,tx,False,1,lon1,lat1,False)

end 
