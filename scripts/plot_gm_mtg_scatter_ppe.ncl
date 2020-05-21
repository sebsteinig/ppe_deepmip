deepmip_dir = "/export/acrc/DeepMIP_Model_Output_shared"
ppe_dir     = "/home/bridge/wb19586/spinup_validation/data/fields"
plot_dir    = "/home/bridge/wb19586/analysis/ppe_deepmip/plots"
grid_ref    = "/home/bridge/wb19586/analysis/ppe_deepmip/tmp/grids/HadCM3BL_M2.1aN-tos.nc"
grid_ref_short    = "HadCM3BL"
;grid_ref    = "r360x180"
;grid_ref_short    = "r360x180"

models_deepmip    = (/"HadCM3BL_M2.1aN","HadCM3B_M2.1aN","GFDL_CM2.1","CESM1.2_CAM5","COSMOS-landveg_r2413","MIROC4m","IPSLCM5A2","NorESM1_F"/)
deepmip_short     = (/"HadCM3BL","HadCM3B","GFDL","CESM","COSMOS","MIROC4m","IPSLCM5A2","NorESM1"/)
versions_deepmip  = (/"v1.0","v1.0","v1.0","v1.0","v1.0","v1.0","v1.0","v1.0"/)
exps_deepmip      = (/"piControl","1x","1.5x","2x","3x","4x","6x","9x"/)
vars_deepmip      = (/"tas","tos","tsurf_gradient","sst_gradient","rmse_tsurf","rmse_sst","pr","pr","psl","TOA_imbalance"/)

models_ppe        = (/"a","b","c","d","e","f","g","h","i","j","k","l","m","n","o","p","q","r","s","t","u","v","w"/)
exps_ppe          = (/"xott","xoud","abc","abc","xoov","xoug","xouh","abc"/)
varnames_ppe      = (/"pd","pg","gradient","gradient","rmse","rmse","pd","pd","pd","TOA"/)
vars_ppe          = (/"temp_mm_1_5m","temp_ym_dpth","tsurf_gradient","sst_gradient","rmse_tsurf","rmse_sst","precip_mm_srf","precip_mm_srf","p_mm_msl","TOA_imbalance"/)
vars_ppe_fields   = (/"temp_mm_1_5m","temp_ym_dpth","tsurf_gradient","sst_gradient","rmse_tsurf","rmse_sst","precip_mm_srf","precip_mm_srf","p_mm_msl","TOA_imbalance"/)

vars_era5         = (/"t2m","sst","t2m","sst","t2m","sst","tp","tp","msl"/)

calc_deepmip      = 0
calc_ppe          = 0

;--- load ERA5 fields
  fi = addfile(plot_dir+"/../tmp/observations/ERA5_1980-1999.clim.HadCM3BL.nc","r")
  fields_era5  =  new((/dimsizes(vars_era5)+1,12,73,96/),"float")
  means_era5   =  new((/dimsizes(vars_era5)+1/),"float")
  do vv=0,dimsizes(vars_era5)-1
    fields_era5(vv,:,:,:) = short2flt(fi->$vars_era5(vv)$)
  end do
  fields_era5(0:5,:,:,:)  = fields_era5(0:5,:,:,:)-273.15   ; temperatures
  fields_era5(6:7,:,:,:)  = fields_era5(6:7,:,:,:)*1000 ; total precipitation
  fields_era5(8,:,:,:)    = fields_era5(8,:,:,:)/100. ; slp to hPa

  fields_era5_ym          = dim_avg_n_Wrap(fields_era5,1)
  fields_era5_ym_zm       = dim_avg_n_Wrap(fields_era5_ym,2)

  lat      := tofloat(fields_era5_ym_zm&latitude)
  ;--- calculate means (area weights from sine from latitude edges)
    pi                            = get_pi("float" )
    dlat                         := new(dimsizes(lat),"float") ; weights
    lats_edge                    := new(dimsizes(lat)+1,"float") ; weights
    if (lat(0) .eq. max(lat)) then
      lats_edge(0)                  = pi/2.0 
      lats_edge(dimsizes(lat))      = -1.*pi/2.0
    else if (lat(0) .eq. min(lat)) then
      lats_edge(0)                  = -1.*pi/2.0 
      lats_edge(dimsizes(lat))      = pi/2.0
    end if
    end if             
    lats_edge(1:dimsizes(lat)-1)  = 2.*pi*(lat(0:dimsizes(lat)-2) + lat(1:dimsizes(lat)-1))/(2.0*360.)
    
    do vv=0,dimsizes(vars_era5)-1
      dlat                          := sin(lats_edge(1:dimsizes(lat))) - sin(lats_edge(0:dimsizes(lat)-1))
      dlat                           = where(ismissing(fields_era5_ym_zm(vv,:)),dlat@_FillValue,dlat)
      if (vv .eq. 2 .or. vv .eq. 3) then
        continue
      end if
      ;--- global means 
        means_era5(vv)               = sum(fields_era5_ym_zm(vv,:)*dlat)/sum(dlat)
      ;--- meridional temperature gradient
        if (vars_era5(vv) .eq. "t2m" .or. vars_era5(vv) .eq. "sst") then
          dlat_low        = where(abs(lat) .lt. 30., dlat,0)
          dlat_high       = where(abs(lat) .gt. 60., dlat,0)
          means_era5(vv+2)  =  wgt_areaave(fields_era5_ym(vv,:,:),dlat_low,1,0) - wgt_areaave(fields_era5_ym(vv,:,:),dlat_high,1,0) ; weighted MTG
        end if
    end do

;--- load DeepMIP TAS+TOS and calculate means
  if (calc_deepmip .eq. 1) then
    means_deepmip     = new((/dimsizes(models_deepmip),dimsizes(exps_deepmip),dimsizes(vars_deepmip)/),"float")
    do nn=0,dimsizes(models_deepmip)-1
      do ii=0,dimsizes(exps_deepmip)-1
        do vv=0,dimsizes(vars_deepmip)-1
        
          if (vars_deepmip(vv) .eq. "TOA_imbalance") then
            if (models_deepmip(nn) .eq. "COSMOS-landveg_r2413" .or. models_deepmip(nn) .eq. "CESM1.2_CAM5") then ; wrong variables
              continue
            end if
            f_rsdt := systemfunc("ls -r "+deepmip_dir+"/*/"+models_deepmip(nn)+"/*"+exps_deepmip(ii)+"*/"+versions_deepmip(nn)+"/*rsdt-"+versions_deepmip(nn)+".mean.nc* 2> /dev/null")
            f_rsut := systemfunc("ls -r "+deepmip_dir+"/*/"+models_deepmip(nn)+"/*"+exps_deepmip(ii)+"*/"+versions_deepmip(nn)+"/*rsut-"+versions_deepmip(nn)+".mean.nc* 2> /dev/null")
            f_rlut := systemfunc("ls -r "+deepmip_dir+"/*/"+models_deepmip(nn)+"/*"+exps_deepmip(ii)+"*/"+versions_deepmip(nn)+"/*rlut-"+versions_deepmip(nn)+".mean.nc* 2> /dev/null")
            if (.not. all(ismissing(f_rsdt))) then ; file(s) found
              remap_file  = plot_dir+"/../tmp/"+models_deepmip(nn)+"_"+exps_deepmip(ii)+"_TOA_imbalance_"+grid_ref_short+".nc"
              if (fileexists(remap_file) .eq. False) then
                calc_TOA_1  = systemfunc("cdo -r sub "+f_rsdt+" "+f_rsut+" "+plot_dir+"/../tmp/TOA.tmp1.nc")
                calc_TOA_2  = systemfunc("cdo -r sub "+plot_dir+"/../tmp/TOA.tmp1.nc"+" "+f_rlut+" "+plot_dir+"/../tmp/TOA.tmp2.nc")
                remap       = systemfunc("cdo -remapnn,"+grid_ref+" "+plot_dir+"/../tmp/TOA.tmp2.nc"+" "+remap_file)
              end if
              fi          = addfile(remap_file,"r")
              data        = fi->rsdt   
            else
              continue
            end if
          else
          ;--- check whether file for model+exp combination is available
            f := systemfunc("ls -r "+deepmip_dir+"/*/"+models_deepmip(nn)+"/*"+exps_deepmip(ii)+"*/"+versions_deepmip(nn)+"/*"+vars_deepmip(vv)+"-"+versions_deepmip(nn)+".mean.nc* 2> /dev/null")
            if (.not. all(ismissing(f))) then ; file(s) found
              remap_file  = plot_dir+"/../tmp/"+models_deepmip(nn)+"_"+exps_deepmip(ii)+"_"+vars_deepmip(vv)+"_"+grid_ref_short+".nc"
              if (fileexists(remap_file) .eq. False) then
                remap       = systemfunc("cdo -remapnn,"+grid_ref+" "+f(0)+" "+remap_file)
              end if
              fi          = addfile(remap_file,"r")
              data        = fi->$vars_deepmip(vv)$
              if (models_deepmip(nn) .eq. "MIROC4m" .and. vars_deepmip(vv) .eq. "psl") then
                data      = data*100.
              end if
            else
              continue
            end if
          end if
          
            if (vars_deepmip(vv) .eq. "tas") then
              data      = data-273.15
            else if (vars_deepmip(vv) .eq. "pr") then
              data      = data*86400
            else if (vars_deepmip(vv) .eq. "psl") then
              data      = data/100.
            end if
            end if
            end if
            if (typeof(data) .eq. "double") then
              data_ym     = dim_avg_n_Wrap(dble2flt(data),0)
            else
              data_ym     = dim_avg_n_Wrap(data,0)
            end if
            data_ym_zm   = dim_avg_n_Wrap(rm_single_dims(data_ym),1) ; calculate annual+zonal mean and remove vertical dimension if present
            data_ym_zm!0 = "lat"
            lat      := tofloat(data_ym_zm&lat)
            ;--- calculate means (area weights from sine from latitude edges)
              pi                            = get_pi("float" )
              dlat                         := new(dimsizes(lat),"float") ; weights
              lats_edge                    := new(dimsizes(lat)+1,"float") ; weights
              if (lat(0) .eq. max(lat)) then
                lats_edge(0)                  = pi/2.0 
                lats_edge(dimsizes(lat))      = -1.*pi/2.0
              else if (lat(0) .eq. min(lat)) then
                lats_edge(0)                  = -1.*pi/2.0 
                lats_edge(dimsizes(lat))      = pi/2.0
              end if
              end if             
              lats_edge(1:dimsizes(lat)-1)  = 2.*pi*(lat(0:dimsizes(lat)-2) + lat(1:dimsizes(lat)-1))/(2.0*360.)
              dlat                          = sin(lats_edge(1:dimsizes(lat))) - sin(lats_edge(0:dimsizes(lat)-1))
              dlat                          = where(ismissing(data_ym_zm),dlat@_FillValue,dlat)
              ;--- global means 
                means_deepmip(nn,ii,vv)     = wgt_areaave(data_ym,dlat,1,0)
              ;--- meridional temperature gradient
                if (vars_deepmip(vv) .eq. "tas" .or. vars_deepmip(vv) .eq. "tos") then
                  dlat_low        = where(abs(lat) .lt. 30., dlat,0)
                  dlat_high       = where(abs(lat) .gt. 60., dlat,0)
                  ;means_deepmip(nn,ii,vv+2)     = sum(data_ym_zm*dlat_low)/sum(dlat_low) - sum(data_ym_zm*dlat_high)/sum(dlat_high) ; weighted MTG from zm
                  means_deepmip(nn,ii,vv+2)     = wgt_areaave(data_ym,dlat_low,1,0) - wgt_areaave(data_ym,dlat_high,1,0) ; weighted MTG
                  delete([/dlat_low,dlat_high/])       
                end if
                
                if (vars_deepmip(vv) .eq. "pr" .and. vars_deepmip(vv+1) .eq. "pr") then
                  delete([/data,data_ym_zm/])
                  continue
                end if
             
                if (vars_deepmip(vv) .eq. "tas" .or. vars_deepmip(vv) .eq. "tos" .or. vars_deepmip(vv) .eq. "pr" .or. vars_deepmip(vv) .eq. "uas" .or. vars_deepmip(vv) .eq. "psl") then
                  if (ii .eq. 0) then ; calculate RMSE only for PI
                  ;--- annual mean RMSE vs. ERA5
                    if (vars_deepmip(vv) .eq. "tas" .or. vars_deepmip(vv) .eq. "tos") then
                      index = vv+4
                    else
                      index = vv
                    end if
                    lat_era   = fields_era5_ym&latitude
                    lat_model = data_ym&latitude
                    lon_era   = fields_era5_ym&longitude
                    lon_model = data_ym&longitude
                    if (all(lat_era .eq. lat_model) .and. all(lon_era .eq. lon_model)) then
                      means_deepmip(nn,ii,index)    = wgt_arearmse(data_ym,fields_era5_ym(vv,:,:),dlat,1.0,0)
                    else 
                      print("grids of reanalysis and model are different - can not calculate RMSE. EXIT!")
                      exit
                    end if
                  end if 
                end if
            delete([/data,data_ym_zm/])         
        end do
      end do
    end do

    ;---Open a new NetCDF file to write to
      fout_name = plot_dir+"/../tmp/DeepMIP_processed_means."+grid_ref_short+".nc"
      system("rm -f " + fout_name)
      fout = addfile(fout_name,"c")
      means_deepmip@long_name = "processed DeepMIP metrics"
      means_deepmip!0         = "model"
      means_deepmip!1         = "experiment"
      means_deepmip!2         = "variable"
      fout->means_deepmip = means_deepmip
  else
    f             := addfile(plot_dir+"/../tmp/DeepMIP_processed_means."+grid_ref_short+".nc","r")
    means_deepmip  = f->means_deepmip
  end if

;--- load HadCM3 PPE runs and calculate means
  if (calc_ppe .eq. 1) then
    means_ppe     = new((/dimsizes(models_ppe),dimsizes(exps_ppe),dimsizes(vars_ppe)/),"float")
    do nn=0,dimsizes(models_ppe)-1
      print("processing "+models_ppe(nn))
      do ii=0,dimsizes(exps_ppe)-1
        do vv=0,dimsizes(vars_ppe)-1
          ;--- check whether file for model+exp combination is available
          f := systemfunc("ls -r "+ppe_dir+"/"+exps_ppe(ii)+models_ppe(nn)+"*"+varnames_ppe(vv)+"*.mean.zm.nc* 2> /dev/null")
          if (.not. all(ismissing(f))) then ; file(s) found
            fi_zm         = addfile(f,"r")
            data_ym_zm    = fi_zm->$vars_ppe(vv)$(0,0,:,0)
            data_ym_zm!0  = "lat"
            lat          := tofloat(data_ym_zm&lat)

            f        = systemfunc("ls -r "+ppe_dir+"/"+exps_ppe(ii)+models_ppe(nn)+"*"+varnames_ppe(vv)+"*.mean.nc* 2> /dev/null")
            fi      := addfile(f,"r")
            data     = fi->$vars_ppe_fields(vv)$(:,0,:,:)
            data_ym  = dim_avg_n_Wrap(data,0)
            if (vars_ppe_fields(vv) .eq. "temp_mm_1_5m") then
              data_ym      = data_ym-273.15
            else if (vars_ppe_fields(vv) .eq. "precip_mm_srf") then
              data_ym      = data_ym*86400
            else if (vars_ppe_fields(vv) .eq. "p_mm_msl") then
              data_ym      = data_ym/100.
            end if
            end if
            end if
    
            ;--- calculate means (area weights from sine from latitude edges)
              pi                            = get_pi("float" )
              dlat                         := new(dimsizes(lat),"float") ; weights
              lats_edge                    := new(dimsizes(lat)+1,"float") ; weights
              if (lat(0) .eq. max(lat)) then
                lats_edge(0)                  = pi/2.0 
                lats_edge(dimsizes(lat))      = -1.*pi/2.0
              else if (lat(0) .eq. min(lat)) then
                lats_edge(0)                  = -1.*pi/2.0 
                lats_edge(dimsizes(lat))      = pi/2.0
              end if
              end if             
              lats_edge(1:dimsizes(lat)-1)  = 2.*pi*(lat(0:dimsizes(lat)-2) + lat(1:dimsizes(lat)-1))/(2.0*360.)
              dlat                          = sin(lats_edge(1:dimsizes(lat))) - sin(lats_edge(0:dimsizes(lat)-1))
              dlat                          = where(ismissing(data_ym_zm),dlat@_FillValue,dlat)
              ;--- global means 
                means_ppe(nn,ii,vv)     = wgt_areaave(data_ym,dlat,1,0)
                  
                if (vars_ppe(vv) .eq. "temp_mm_1_5m") ; calculate TOA imbalance
                  SW_down     = fi_zm->downSol_mm_TOA(0,0,:,0)
                  SW_up       = fi_zm->upSol_mm_s3_TOA(0,0,:,0)
                  LW_up       = fi_zm->olr_mm_s3_TOA(0,0,:,0)
                  TOA_balance = SW_down-SW_up-LW_up
                  means_ppe(nn,ii,9)      = sum(TOA_balance*dlat)/sum(dlat)
                end if

                if (vars_ppe_fields(vv) .eq. "precip_mm_srf" .and. vars_ppe_fields(vv+1) .eq. "precip_mm_srf") then
                  continue
                end if

              ;--- SST gradient
                if (vars_ppe(vv) .eq. "temp_ym_dpth" .or. vars_ppe(vv) .eq. "temp_mm_1_5m") then
                  dlat_low        = where(abs(lat) .lt. 30., dlat,0)
                  dlat_high       = where(abs(lat) .gt. 60., dlat,0)
                  ;means_ppe(nn,ii,vv+2)     = sum(data_ym_zm*dlat_low)/sum(dlat_low) - sum(data_ym_zm*dlat_high)/sum(dlat_high) ; weighted MTG from zm
                  means_ppe(nn,ii,vv+2)     =  wgt_areaave(data_ym,dlat_low,1,0) - wgt_areaave(data_ym,dlat_high,1,0); weighted MTG
                  delete([/dlat_low,dlat_high/])      
                end if
              ;--- annual mean RMSE vs. ERA5
              if (vars_ppe_fields(vv) .eq. "temp_ym_dpth" .or. vars_ppe_fields(vv) .eq. "temp_mm_1_5m" .or. vars_ppe_fields(vv) .eq. "precip_mm_srf" .or. vars_ppe_fields(vv) .eq. "u_mm_10m" .or. vars_ppe_fields(vv) .eq. "p_mm_msl") then
                  if (ii .eq. 0) then ; calculate RMSE only for PI
                    if (vars_ppe_fields(vv) .eq. "temp_mm_1_5m" .or. vars_ppe_fields(vv) .eq. "temp_ym_dpth") then
                      index = vv+4
                    else
                      index = vv
                    end if
                
                    lat_era   = fields_era5_ym&latitude
                    lat_model = data_ym&latitude
                    lon_era   = fields_era5_ym&longitude
                    lon_model = data_ym&longitude
                    if (all(lat_era .eq. lat_model) .and. all(lon_era .eq. lon_model)) then
                      means_ppe(nn,ii,index)    = wgt_arearmse(data_ym,fields_era5_ym(vv,:,:),dlat,1.0,0)
                    else if (all(lat_era .eq. lat_model(::-1)) .and. all(lon_era .eq. lon_model)) then ; invert latitude of HadCM3 atmospheric data
                      means_ppe(nn,ii,index)    = wgt_arearmse(data_ym(::-1,:),fields_era5_ym(vv,:,:),dlat,1.0,0)
                    else 
                      print("grids of reanalysis and model are different - can not calculate RMSE. EXIT!")
                      exit
                    end if
                    end if
                  end if
                end if
            delete([/data_ym_zm,data,data_ym/])         
          end if
        end do
      end do
    end do

    ;---Open a new NetCDF file to write to
      fout_name = plot_dir+"/../tmp/PPE_processed_means.HadCM3BL.nc"
      system("rm -f " + fout_name)
      fout = addfile(fout_name,"c")
      means_ppe@long_name = "processed HadCM3B PPE metrics"
      means_ppe!0         = "model"
      means_ppe!1         = "experiment"
      means_ppe!2         = "variable"
      fout->means_ppe     = means_ppe
  else
    f             := addfile(plot_dir+"/../tmp/PPE_processed_means.HadCM3BL.nc","r")
    means_ppe      = f->means_ppe
  end if

;-- mask PPE means with TOA imbalance
  means_ppe = where(conform_dims(dimsizes(means_ppe),abs(means_ppe(:,:,9)),(/0,1/)) .gt. 0.5,means_ppe@_FillValue,means_ppe)

;--- create mask for PI constraints
  means_ppe_mask  = means_ppe
  means_ppe_mask  = 0

  min_deepmip_all   = dim_avg_n(means_deepmip(2:7,0,0:8),0)-2*dim_stddev_n(means_deepmip(2:7,0,0:8),0)
  max_deepmip_all   = dim_avg_n(means_deepmip(2:7,0,0:8),0)+2*dim_stddev_n(means_deepmip(2:7,0,0:8),0)

  do nn=0,dimsizes(models_ppe)-1
    flag=1
    do vv=0,8
      if (vv .ne. 8) then
        if (means_ppe(nn,0,vv) .lt. min_deepmip_all(vv) .or. means_ppe(nn,0,vv) .gt. max_deepmip_all(vv)) then
          flag=0
        end if
      end if
    end do
    means_ppe_mask(nn,:,:) = flag
  end do

print(models_ppe+": "+means_ppe_mask(:,0,0))
;--- scatter plots
  wks = gsn_open_wks("pdf",plot_dir+"/gm_mtg_scatter_ppe")

  res_blank                 = True
  res_blank@gsnDraw         = False
  res_blank@gsnFrame        = False
  res_blank@tmXBMode        = "Explicit"
  res_blank@vpWidthF        = 0.6
  res_blank@vpHeightF       = 0.4

  res                       = res_blank
  res@xyMarkLineMode        = "Markers"
  res@gsnFrame              = False
  
  resP                = True
  resP@gsnMaximize    = True
  resP@gsnPanelYWhiteSpacePercent = 5.
  resP@gsnPanelXWhiteSpacePercent = 5.

  plres                 = True
  plres@gsMarkerSizeF   = 0.012
  
  ploutres                    = True
  ploutres@gsMarkerSizeF      = 0.007
  ploutres@gsMarkerThicknessF = 1.
  ploutres@gsMarkerIndex      = 4

  pl_colors             = (/"gray","black","black","purple","blue","cyan","green","orange"/)
  pl_styles_deepmip     = (/4,16,5,6,9,7,8,3/)
  pl_styles_ppe_f       = (/12/)
  pl_styles_ppe         = (/16/)

  pgres                 = True
  pgres@gsFillColor     = "transparent"      ;-- fill color
  pgres@gsEdgesOn       = True
  pgres@gsEdgeColor     = "red"         ;-- set opacity of polygon
  pgres@gsEdgeThicknessF = 2.0
  
  plot          = new(4,graphic)
  poly          = new(2000,graphic)
  poly2         = new(2000,graphic)
  poly2_out     = new(2000,graphic)

  count         = 0
  
  co2_levels    = (/1,1,1.5,2,3,4,6,9/)
  
    ;--- TEMP2 vs. CO2
    res_blank@trXMinF         = log10(0.7)
    res_blank@trXMaxF         = log10(20)
    res_blank@trYMinF         = 10
    res_blank@trYMaxF         = 37
    res_blank@tmXBValues     = log10((/1,2,4,8,16/))
    res_blank@tmXBLabels     = (/"1x","2x","4x","8x","16x"/)
    res_blank@tiXAxisString  = "CO2 relative to preindustrial"
    res_blank@tiYAxisString  = "global mean surface temperature"
    res_blank@tiXAxisFontHeightF  = 0.016
    res_blank@tiYAxisFontHeightF  = 0.016
    res_blank@tmXBLabelFontHeightF  = 0.015
    res_blank@tmYLLabelFontHeightF  = 0.015

    plot(0)       = gsn_csm_blank_plot(wks,res_blank)
    
    do nn=0,dimsizes(models_deepmip)-1
      plres@gsMarkerSizeF   = 0.012
      plres@gsMarkerThicknessF = 1.5
      do ll=0,dimsizes(co2_levels)-1
        plres@gsMarkerIndex = pl_styles_deepmip(nn)
        plres@gsMarkerColor = pl_colors(ll)
        poly(count) = gsn_add_polymarker(wks,plot(0),log10(co2_levels(ll)-0.05*co2_levels(ll)),means_deepmip(nn,ll,0),plres)
        count       = count+1
        plres@gsMarkerIndex   = pl_styles_ppe_f
        poly(count) = gsn_add_polymarker(wks,plot(0),log10(co2_levels(ll)+0.05*co2_levels(ll)),means_ppe(5,ll,0),plres)
        count       = count+1
      end do
    end do   
    
    proxy_area_a  = gsn_add_polygon(wks,plot(0),log10((/3.4,9,9,3.4,3.4/)),(/26,26,32,32,26/),pgres)

    plot(2)       = gsn_csm_blank_plot(wks,res_blank)

    plres@gsMarkerOpacityF = 1.0
    do nn=0,dimsizes(models_deepmip)-1
      plres@gsMarkerSizeF   = 0.012
      plres@gsMarkerThicknessF = 1.5
      plres@gsMarkerIndex = pl_styles_deepmip(nn)
      do ll=0,dimsizes(co2_levels)-1
        plres@gsMarkerColor = pl_colors(ll)
        poly2(count) = gsn_add_polymarker(wks,plot(2),log10(co2_levels(ll)-0.05*co2_levels(ll)),means_deepmip(nn,ll,0),plres)
        count       = count+1
      end do
    end do   
    plres@gsMarkerOpacityF = 1.
    do nn=0,dimsizes(models_ppe)-1
      if (models_ppe(nn) .eq. "f") then
        plres@gsMarkerIndex   = pl_styles_ppe_f
        plres@gsMarkerSizeF   = 0.012
        plres@gsMarkerThicknessF = 1.5
      else
        plres@gsMarkerSizeF   = 0.007
        plres@gsMarkerThicknessF = 1.
        plres@gsMarkerIndex = pl_styles_ppe
      end if
      do ll=0,dimsizes(co2_levels)-1
        plres@gsMarkerColor = pl_colors(ll)
        
        if (means_ppe_mask(nn,ll,0) .eq. 1) then
          poly2(count) = gsn_add_polymarker(wks,plot(2),log10(co2_levels(ll)+0.05*co2_levels(ll)),means_ppe(nn,ll,0),plres)
        end if
        if (models_ppe(nn) .ne. "f") then
          if (means_ppe_mask(nn,ll,0) .eq. 1) then
            ploutres@gsMarkerColor  = "black"
          else
            ploutres@gsMarkerColor  = pl_colors(ll)
          end if
          poly2_out(count) = gsn_add_polymarker(wks,plot(2),log10(co2_levels(ll)+0.05*co2_levels(ll)),means_ppe(nn,ll,0),ploutres)
        end if
        count       = count+1
      end do
    end do
    
    proxy_area_a2  = gsn_add_polygon(wks,plot(2),log10((/3.4,9,9,3.4,3.4/)),(/26,26,32,32,26/),pgres)

    ;--- SST vs. SST MTG
    res_blank@trXMinF         = 15
    res_blank@trXMaxF         = 40
    res_blank@trYMinF         = 16
    res_blank@trYMaxF        := 28.5
    res_blank@tmXBValues     := ispan(15,40,5)
    res_blank@tmXBLabels     := (/"15","20","25","30","35","40"/)
    res_blank@tmXBMinorValues:= ispan(15,40,1)

    res_blank@tiXAxisString  = "global mean SST"
    res_blank@tiYAxisString  = "meridional SST gradient"

    plot(1)       = gsn_csm_blank_plot(wks,res_blank)
    
    do nn=0,dimsizes(models_deepmip)-1
      plres@gsMarkerSizeF   = 0.012
      plres@gsMarkerThicknessF = 1.5
      do ll=0,dimsizes(co2_levels)-1
        plres@gsMarkerIndex = pl_styles_deepmip(nn)
        plres@gsMarkerColor = pl_colors(ll)
        poly(count) = gsn_add_polymarker(wks,plot(1),means_deepmip(nn,ll,1),means_deepmip(nn,ll,3),plres)
        count       = count+1
        plres@gsMarkerIndex   = pl_styles_ppe_f
        poly(count) = gsn_add_polymarker(wks,plot(1),means_ppe(5,ll,1),means_ppe(5,ll,3),plres)
        count       = count+1        
      end do
    end do

    proxy_area_b  = gsn_add_polygon(wks,plot(1),(/26,32,32,26,26/),(/20,20,22,22,20/),pgres)

    plot(3)       = gsn_csm_blank_plot(wks,res_blank)
    
    plres@gsMarkerOpacityF = 0.3
    do nn=0,dimsizes(models_deepmip)-1
      plres@gsMarkerSizeF   = 0.012
      plres@gsMarkerThicknessF = 2.
      plres@gsMarkerIndex = pl_styles_deepmip(nn)
      do ll=0,dimsizes(co2_levels)-1
        plres@gsMarkerColor = pl_colors(ll)
        poly2(count) = gsn_add_polymarker(wks,plot(3),means_deepmip(nn,ll,1),means_deepmip(nn,ll,3),plres)
        count       = count+1
      end do
    end do
    plres@gsMarkerOpacityF = 1.0
    do nn=0,dimsizes(models_ppe)-1
      if (models_ppe(nn) .eq. "f") then
        plres@gsMarkerIndex = pl_styles_ppe_f
        plres@gsMarkerSizeF   = 0.012
        plres@gsMarkerThicknessF = 2.
      else
        plres@gsMarkerSizeF   = 0.007
        plres@gsMarkerThicknessF = 1.
        plres@gsMarkerIndex = pl_styles_ppe
      end if
      do ll=0,dimsizes(co2_levels)-1
        plres@gsMarkerColor = pl_colors(ll)
        if (means_ppe_mask(nn,ll,0) .eq. 1) then
          poly2(count) = gsn_add_polymarker(wks,plot(3),means_ppe(nn,ll,1),means_ppe(nn,ll,3),plres)
        end if
        if (models_ppe(nn) .ne. "f") then
          if (means_ppe_mask(nn,ll,0) .eq. 1) then
            ploutres@gsMarkerColor  = "black"
          else
            ploutres@gsMarkerColor  = pl_colors(ll)
          end if
          poly2_out(count) = gsn_add_polymarker(wks,plot(3),means_ppe(nn,ll,1),means_ppe(nn,ll,3),ploutres)
        end if
        count       = count+1
      end do
    end do

    proxy_area_b2  = gsn_add_polygon(wks,plot(3),(/26,32,32,26,26/),(/20,20,22,22,20/),pgres)

   names_legend = array_append_record(models_deepmip,"HadCM3B_M2.1aN_PPE_f",0)
   pl_styles_legend = array_append_record(pl_styles_deepmip,pl_styles_ppe_f,0)
   
   legend = create "Legend" legendClass wks 
   "vpXF" : 0.305
   "vpYF" : .665
   "lgAutoManage" : False
   "vpWidthF" : 0.08
   "vpHeightF" : 0.12
   "lgLineThicknessF" : 8.0
   "lgLineDashSegLenF" : 0.08
   "lgPerimOn" : True
   "lgLabelJust" : "CenterLeft"
   "lgItemCount" : dimsizes(names_legend)
   "lgItemOrder" : ispan(dimsizes(names_legend)-1,0,1)
   "lgLabelStrings" : names_legend
   "lgItemType" : "Markers"
   "lgLabelsOn" : True 
   "lgLineLabelsOn" : False 
   "lgLabelFontHeightF"  : 0.006
   "lgLabelPosition" : "Right"
   "lgLabelOffsetF" : 0.1
   "lgMonoDashIndex" : False
   "lgMarkerIndexes" : pl_styles_legend
   "lgMonoMarkerColor" : True ; Only one dash pattern - solid
   "lgMarkerColor" : "black"
   "lgMarkerSizeF" : 0.008
   "lgMarkerThicknessF" : 2
   end create 

   draw(legend)

   names_legend2 = (/"TOA < 0.5","TOA < 0.5 + 'good' PI"/)
   pl_styles_legend2 = (/4,pl_styles_ppe/)
   
   legend2 = create "Legend" legendClass wks 
   "vpXF" : 0.845
   "vpYF" : .445
   "lgAutoManage" : False
   "vpWidthF" : 0.04
   "vpHeightF" : 0.04
   "lgLineThicknessF" : 8.0
   "lgLineDashSegLenF" : 0.08
   "lgPerimOn" : True
   "lgLabelJust" : "CenterLeft"
   "lgItemCount" : dimsizes(names_legend2)
   "lgItemOrder" : ispan(dimsizes(names_legend2)-1,0,1)
   "lgLabelStrings" : names_legend2
   "lgItemType" : "Markers"
   "lgLabelsOn" : True 
   "lgLineLabelsOn" : False 
   "lgLabelFontHeightF"  : 0.007
   "lgLabelPosition" : "Right"
   "lgLabelOffsetF" : 0.1
   "lgMonoDashIndex" : False
   "lgMarkerIndexes" : pl_styles_legend2
   "lgMonoMarkerColor" : True ; Only one dash pattern - solid
   "lgMarkerColor" : "black"
   "lgMarkerSizeF" : 0.008
   "lgMarkerThicknessF" : 1
   end create 

   draw(legend2)

   names_legend3 = (/"piControl","DeepMIP 1x/1.5x","DeepMIP 2x","DeepMIP 3x","DeepMIP 4x","DeepMIP 6x","DeepMIP 9x"/)
   pl_styles_legend3 = (/pl_styles_ppe,pl_styles_ppe,pl_styles_ppe,pl_styles_ppe,pl_styles_ppe,pl_styles_ppe,pl_styles_ppe/)
   
   legend3 = create "Legend" legendClass wks 
   "vpXF" : 0.08
   "vpYF" : .795
   "lgAutoManage" : False
   "vpWidthF" : 0.03
   "vpHeightF" : 0.1
   "lgLineThicknessF" : 8.0
   "lgLineDashSegLenF" : 0.08
   "lgPerimOn" : False
   "lgLabelJust" : "CenterLeft"
   "lgItemCount" : dimsizes(names_legend3)
   "lgItemOrder" : ispan(dimsizes(names_legend3)-1,0,1)
   "lgLabelStrings" : names_legend3
   "lgItemType" : "Markers"
   "lgLabelsOn" : True 
   "lgLineLabelsOn" : False 
   "lgLabelFontHeightF"  : 0.007
   "lgLabelPosition" : "Right"
   "lgLabelOffsetF" : 0.1
   "lgMonoDashIndex" : False
   "lgMarkerIndexes" : pl_styles_legend3
   "lgMonoMarkerColor" : False ; Only one dash pattern - solid
   "lgMarkerColors" : pl_colors((/0,1,3,4,5,6,7/))
   "lgMarkerSizeF" : 0.008
   "lgMarkerThicknessF" : 1
   end create 

   draw(legend3)

  resP@txString               = "HadCM3B PPE piControl+DeepMIP simulations vs. DeepMIP database vs. EECO proxies (red square) " 
  resP@txFontHeightF          = 0.015
  gsn_panel(wks,plot,(/2,2/),resP)


  wks = gsn_open_wks("pdf",plot_dir+"/PI_metric_comparison")

  idx_list        = (/0,1,2,3,4,5,6,7,8,9/)
  vars_titles     = (/"global mean T2m","global mean SST","meridional gradient T2m","polar amplification","RMSE T2m","RMSE SST","global mean Precip","RMSE Precip","RMSE sea-level pressure","TOA imbalance"/)
  ymins           = (/10,15,  34,23.5,0,0,2.5,0.0 ,0,-.6/)
  ymaxs           = (/18,21.5,47,30,4,4,3.4,2.55,8, 1.2/)

  plot                  := new(dimsizes(idx_list),graphic)
  deepmip_poly          := new(dimsizes(idx_list),graphic)
  era5_poly             := new(dimsizes(idx_list),graphic)
  deepmip_lines         := new(dimsizes(idx_list)*dimsizes(models_deepmip),graphic)
  deepmip_anno          := new(dimsizes(idx_list)*dimsizes(models_deepmip),graphic)

  res_bar                 = True
  res_bar@gsnXYBarChart   = True            ; Create bar plot
  res_bar@gsnDraw         = False
  res_bar@gsnFrame        = False
  res_bar@trXMinF         = -1.
  res_bar@trXMaxF         = 25.
  res_bar@tmXBMode        = "Explicit"
  res_bar@tmXBValues      = ispan(1,23,1)
  res_bar@tmXBLabels      = models_ppe
  res_bar@tmXBLabelJust   = "CenterCenter"
  res_bar@tiYAxisString   = ""
  res_bar@gsnYRefLine     = 0
  colors_bar              = models_ppe
  res_bar@xyCurveDrawOrder  = "PostDraw"
  res_bar@tmXTOn            = False
  res_bar@tmYRLabelsOn    = True  
  res_bar@tmYUseLeft = False
  res_bar@tmYRMode        = "Explicit"  

  gres                 = True
  gres@gsFillColor     = "black"
  gres@gsFillOpacityF  = 0.15
  gres@tfPolyDrawOrder = "Draw"  ; this can be used for polylines, polymarkers, or polygons

  txres                = True
  txres@txFontHeightF  = 0.01
  
  pline_res                   = True

  do nn=0,dimsizes(idx_list)-1
    res_bar@trYMinF       := ymins(nn)
    res_bar@trYMaxF       := ymaxs(nn)
    min_deepmip            = avg(means_deepmip(2:7,0,idx_list(nn)))-2*stddev(means_deepmip(2:7,0,idx_list(nn)))
    max_deepmip            = avg(means_deepmip(2:7,0,idx_list(nn)))+2*stddev(means_deepmip(2:7,0,idx_list(nn)))
    do ii=0,22
      if (means_ppe(ii,0,idx_list(nn)) .lt. min_deepmip .or. means_ppe(ii,0,idx_list(nn)) .gt. max_deepmip) then
        colors_bar(ii)  = "beige"
      else
        if (models_ppe(ii) .eq. "a") then
          colors_bar(ii)  = "gray50"
        else if (models_ppe(ii) .eq. "f") then
          colors_bar(ii)  = "firebrick1"
        else
          colors_bar(ii)  = "burlywood1"
        end if
        end if
      end if
    end do
    if (vars_titles(nn) .eq. "TOA imbalance") then
      delete(res_bar@gsnXYBarChartColors)
      res_bar@gsnXYBarChartColors2 = colors_bar
    else
      res_bar@gsnXYBarChartColors = colors_bar
    end if
    
    res_bar@tmYRValues = means_deepmip(:,0,idx_list(nn))
    res_bar@tmYRLabels = deepmip_short
    res_bar@tmYRLabelFontHeightF = 0.01

    plot(nn)               = gsn_csm_xy(wks,ispan(1,23,1),means_ppe(:,0,idx_list(nn)),res_bar)
    deepmip_poly(nn)       = gsn_add_polygon(wks,plot(nn),(/-10,30,30,-10,-10/),(/min_deepmip,min_deepmip,max_deepmip,max_deepmip,min_deepmip/),gres)
    
    pline_res@gsLineThicknessF  = 1
    pline_res@gsLineColor  = "gray30"
    do dd=0,dimsizes(models_deepmip)-1
      if (.not. ismissing(means_deepmip(dd,0,idx_list(nn)))) then
        deepmip_lines(nn*dimsizes(models_deepmip)+dd)  = gsn_add_polyline(wks,plot(nn),(/-10,30/),(/means_deepmip(dd,0,idx_list(nn)),means_deepmip(dd,0,idx_list(nn))/),pline_res)
      end if
    end do
    if (nn .le. 3 .or. nn .eq. 6) then
      pline_res@gsLineThicknessF  = 3
      pline_res@gsLineColor       = "black"
      era5_poly(nn)          = gsn_add_polyline(wks,plot(nn),(/-10,30/),(/means_era5(idx_list(nn)),means_era5(idx_list(nn))/),pline_res)
    end if
  end do
  
  resP@amJust                 = "TopLeft"
  resP@gsnPanelFigureStrings  = vars_titles
  resP@gsnPanelFigureStringsFontHeightF = 0.01
  resP@gsnPanelXWhiteSpacePercent = 0
  resP@gsnPanelYWhiteSpacePercent = 0
  resP@txString               = "HadCM3B PPE piControl performance vs. DeepMIP database~C~shading: DeepMIP mean +/- 2 sigma; black line and RMSE ref: ERA5 1980-1999 mean " 
  resP@txFontHeightF          = 0.015
  
  gsn_panel(wks,plot,(/3,4/),resP)
