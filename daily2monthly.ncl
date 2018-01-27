load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
begin
;***********************************************
; get variable names from nc files
;***********************************************
   f       = addfile("SESA_prctntx_1961_2000_1deg.nc","r")
   system("rm SESA_prctntx_1961_2000_1deg_monthly.nc") ; remove any pre-existing file
   ncdf_out= addfile("SESA_prctntx_1961_2000_1deg_monthly.nc" ,"c")
   
   pr	   = f->prec
   tn	   = f->tmi
   tx      = f->tma

   
   pr_mon  = calculate_monthly_values(pr, "ave", 0, False)
   tn_mon  = calculate_monthly_values(tn, "ave", 0, False)
   tx_mon  = calculate_monthly_values(tx, "ave", 0, False)

   copy_VarAtts(pr(0,:,:), pr_mon(0,:,:))
   pr_mon@units = "mm/day" 

   filedimdef(ncdf_out,"time",-1,True) 
   time=ispan(0,dimsizes(pr_mon(:,0,0))-1,1)
   filevardef(ncdf_out, "time" ,typeof(pr&time),getvardims(pr&time))
   ;copy_VarAtts(pr&time, time)
   time@units = "months since 1961-01-15 00:00"
   time@long_name="time in months"
   tn_mon@scale_factor = 0.1
   tx_mon@scale_factor = 0.1

   pr_mon@long_name = "monthly precipitation (accumulated)"
   pr_mon&time=time
   tn_mon&time=time
   tx_mon&time=time
   tn_mon@Origin="Barbara Tencer minimum temperature dataset"
   tx_mon@Origin="Barbara Tencer maximum temperature dataset"

   ncdf_out->time  = time
   ncdf_out->prec  = pr_mon
   ncdf_out->tma   = tx_mon
   ncdf_out->tmi   = tn_mon

system("ncdump -h SESA_prctntx_1961_2000_1deg_monthly.nc")
   
end

;--------------------------------- 
;Function begins here!
undef ("calculate_monthly_values") 
function calculate_monthly_values (x:numeric, arith:string, nDim[1]:integer, opt[1]:logical) 
; calculate monthly values [avg, sum, min, max] 
; x: numeric array of 5D or less [eg: time,lev,lat,lon] 
; *must* have time coordinate recognized by ut_calendar 
; if 5D [case,time,lev,lat,lon] 
; arith: "avg" [also, "ave"], "sum","min","max" others may be added later 
; nDim : scalar integer that specifies the 'time' dimension [generally, 0] 
; opt : option ... not used here 
; 
; Sample usage: x(time,lat,lon) where time are n-hrly or daily values. 
; xMonthAvg = calculate_monthly_values(x, "avg", 0, False) 
; xMonthSum = calculate_monthly_values(x, "sum", 0, False) 
; xMonthMin = calculate_monthly_values(x, "min", 0, False) 
; xMonthMax = calculate_monthly_values(x, "max", 0, False) 
; 
local dimx, rankx, utc_date, year, month, day, hour, ntim \ 
    , yrStrt, yrLast, nyrs, NTIM, dAvg, xReturn, xMonth \ 
    , NT, nmo, ii 
begin 
    dimx = dimsizes( x ) 
    rankx = dimsizes( dimx ) 
    if (rankx.gt.5) then 
        print("calculate_monthly_values: rankx="+rankx +" [only 5D or fewer supported]") 
        exit 
    end if 
    utc_date = ut_calendar(x&time, 0) 
    year = floattointeger(utc_date(:,0)) 
    month = floattointeger(utc_date(:,1)) 
    day = floattointeger(utc_date(:,2)) 
    hour = floattointeger(utc_date(:,3)) 
   ;minute = floattointeger(utc_date(:,4)) 
   ;second = utc_date(:,5) 
    if (rankx.le.4) then 
        ntim = dimx(0) 
    else 
        ntim = dimx(1) 
    end if 
    yrStrt = year(0) 
    yrLast = year(ntim-1) 
    nyrs = yrLast-yrStrt+1 
    NTIM = 12*nyrs ; size of monthly files 
    dAvg = dimx 
    if (rankx.le.4) then 
        dAvg(0)= NTIM 
    else 
        dAvg(1)= NTIM 
    end if 
    xReturn= new ( dAvg , typeof(x), getFillValue(x)) 
    xMonth = new ( NTIM , typeof(x&time), "No_FillValue") 
    NT = -1 
    do yr=yrStrt,yrLast 
      do nmo=0,11 
         NT = NT+1 
         if (isvar("ii")) then ; keep this here!! 
             delete(ii) 
         end if 
         ii = ind(yr.eq.year .and. (nmo+1).eq.month) 
         if (.not.ismissing(ii(0))) then 
             xMonth(NT) = (/ x&time(ii(0)) /) 
             if (rankx.eq.1) then 
                 if (arith.eq."avg" .or. arith.eq."ave") then 
                     xReturn(NT) = dim_avg_n(x(ii) , nDim) 
                 end if 
                 if (arith.eq."sum") then 
                     xReturn(NT) = dim_sum_n(x(ii) , nDim) 
                 end if 
                 if (arith.eq."min") then 
                     xReturn(NT) = dim_min_n(x(ii) , nDim) 
                 end if 
                 if (arith.eq."max") then 
                     xReturn(NT) = dim_max_n(x(ii) , nDim) 
                 end if 
             end if 
             if (rankx.eq.2) then 
                 if (arith.eq."avg" .or. arith.eq."ave") then 
                     xReturn(NT,:) = dim_avg_n(x(ii,:) , nDim) 
                 end if 
                 if (arith.eq."sum") then 
                     xReturn(NT,:) = dim_sum_n(x(ii,:) , nDim) 
                 end if 
                 if (arith.eq."min") then 
                     xReturn(NT,:) = dim_min_n(x(ii,:) , nDim) 
                 end if 
                 if (arith.eq."max") then 
                     xReturn(NT,:) = dim_max_n(x(ii,:) , nDim) 
                 end if 
             end if 
             if (rankx.eq.3) then 
                 if (arith.eq."avg" .or. arith.eq."ave") then 
                     xReturn(NT,:,:) = dim_avg_n(x(ii,:,:) , nDim) 
                 end if 
                 if (arith.eq."sum") then 
                     xReturn(NT,:,:) = dim_sum_n(x(ii,:,:) , nDim) 
                 end if 
                 if (arith.eq."min") then 
                     xReturn(NT,:,:) = dim_min_n(x(ii,:,:) , nDim) 
                 end if 
                 if (arith.eq."max") then 
                     xReturn(NT,:,:) = dim_max_n(x(ii,:,:) , nDim) 
                 end if 
             end if 
             if (rankx.eq.4) then 
                 if (arith.eq."avg" .or. arith.eq."ave") then 
                     xReturn(NT,:,:,:) = dim_avg_n(x(ii,:,:,:) , nDim) 
                 end if 
                 if (arith.eq."sum") then 
                     xReturn(NT,:,:,:) = dim_sum_n(x(ii,:,:,:) , nDim) 
                 end if 
                 if (arith.eq."min") then 
                     xReturn(NT,:,:,:) = dim_min_n(x(ii,:,:,:) , nDim) 
                 end if 
                 if (arith.eq."max") then 
                     xReturn(NT,:,:,:) = dim_max_n(x(ii,:,:,:) , nDim) 
                 end if 
             end if 
             if (rankx.eq.5) then ; note the location of time 
                 if (arith.eq."avg" .or. arith.eq."ave") then 
                     xReturn(:,NT,:,:,:) = dim_avg_n(x(:,ii,:,:,:) , nDim) 
                 end if 
                 if (arith.eq."sum") then 
                     xReturn(:,NT,:,:,:) = dim_sum_n(x(:,ii,:,:,:) , nDim) 
                 end if 
                 if (arith.eq."min") then 
                     xReturn(:,NT,:,:,:) = dim_min_n(x(:,ii,:,:,:) , nDim) 
                 end if 
                 if (arith.eq."max") then 
                     xReturn(:,NT,:,:,:) = dim_max_n(x(:,ii,:,:,:) , nDim) 
                 end if 
             end if 
         end if 
         delete(ii) 
      end do ; month 
    end do ; year 
    if (rankx.le.4) then 
        xReturn!0= "time" 
    else 
        xReturn!1= "time" 
    end if 
    xMonth@units = x&time@units 
    xReturn&time = xMonth 
    if (isatt(x,"long_name")) then 
        xReturn@long_name = "Monthly precipitation (accumulated)" ;x@long_name 
    end if 
    if (isatt(x,"units")) then 
        xReturn@units = x@units 
    end if 
    if (rankx.eq.3) then 
        copy_VarCoords(x(0,:,:), xReturn(0,:,:)) 
    end if 
    if (rankx.eq.4) then 
        copy_VarCoords(x(0,0,:,:), xReturn(0,0,:,:)) 
    end if 
    if (rankx.eq.5) then 
        copy_VarCoords(x(:,0,0,:,:), xReturn(:,0,0,:,:)) 
    end if 
    xReturn@operation_tag = "calculate_monthly_values: "+arith 
    return( xReturn ) 
end 