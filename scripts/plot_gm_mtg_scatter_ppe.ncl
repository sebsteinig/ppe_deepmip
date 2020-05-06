deepmip_dir = "/export/acrc/DeepMIP_Model_Output_shared"
ppe_dir     = "/home/bridge/wb19586/spinup_validation/data/fields"
plot_dir    = "/home/bridge/wb19586/analysis/ppe_deepmip/plots"
grid_ref    = "/home/bridge/wb19586/analysis/ppe_deepmip/tmp/grids/HadCM3BL_M2.1aN-tos.nc"
grid_ref_short    = "HadCM3BL"
;grid_ref    = "r360x180"
;grid_ref_short    = "r360x180"

models_deepmip    = (/"HadCM3BL_M2.1aN","HadCM3B_M2.1aN","GFDL_CM2.1","CESM1.2_CAM5","COSMOS-landveg_r2413","MIROC4m"/)
versions_deepmip  = (/"v1.0","v1.0","v1.0","v1.0","v1.0","v1.0"/)
exps_deepmip      = (/"{pi_control,piControl}","1x","2x","3x","4x","6x","9x"/)
vars_deepmip      = (/"tas","tos"/)

models_ppe        = (/"a","b","c","d","e","f","g","h","i","j","k","l","m","n","o","p","q","r","s","t","u","v","w"/)
exps_ppe          = (/"xott","xoud","abc","xoov","abc","xouh","abc"/)
varnames_ppe      = (/"pd","pg"/)
vars_ppe          = (/"temp_mm_1_5m","temp_ym_dpth"/)

calc_deepmip      = 0
calc_ppe          = 0

;--- load DeepMIP TAS+TOS and calculate means
  if (calc_deepmip .eq. 1) then
    means_deepmip     = new((/dimsizes(models_deepmip),dimsizes(exps_deepmip),dimsizes(vars_deepmip)+1/),"float")
    do nn=0,dimsizes(models_deepmip)-1
      do ii=0,dimsizes(exps_deepmip)-1
        do vv=0,dimsizes(vars_deepmip)-1
          ;--- check whether file for model+exp combination is available
          f := systemfunc("ls -r "+deepmip_dir+"/*/"+models_deepmip(nn)+"/*"+exps_deepmip(ii)+"*/"+versions_deepmip(nn)+"/*"+vars_deepmip(vv)+"-"+versions_deepmip(nn)+".mean.nc* 2> /dev/null")
          if (.not. all(ismissing(f))) then ; file(s) found
            if (ii .eq. 0) then
              remap_file  = plot_dir+"/../tmp/"+models_deepmip(nn)+"_piControl_"+vars_deepmip(vv)+"_"+grid_ref_short+".nc"
            else
              remap_file  = plot_dir+"/../tmp/"+models_deepmip(nn)+"_"+exps_deepmip(ii)+"_"+vars_deepmip(vv)+"_"+grid_ref_short+".nc"
            end if
            if (fileexists(remap_file) .eq. False) then
              remap       = systemfunc("cdo -remapnn,"+grid_ref+" "+f(0)+" "+remap_file)
            end if
            fi          = addfile(remap_file,"r")
            data        = fi->$vars_deepmip(vv)$
            data_ym_zm   = dim_avg_n_Wrap(rm_single_dims(dim_avg_n_Wrap(data,0)),1) ; calculate annual+zonal mean and remove vertical dimension if present
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
                means_deepmip(nn,ii,vv)     = sum(data_ym_zm*dlat)/sum(dlat)
              ;--- SST gradient
                if (vars_deepmip(vv) .eq. "tos") then
                  data_ym_zm_low  = where(abs(lat) .lt. 30., data_ym_zm,data_ym_zm@_FillValue)
                  data_ym_zm_high = where(abs(lat) .gt. 60., data_ym_zm,data_ym_zm@_FillValue)
                  dlat_low        = where(abs(lat) .lt. 30., dlat,dlat@_FillValue)
                  dlat_high       = where(abs(lat) .gt. 60., dlat,dlat@_FillValue)

                  means_deepmip(nn,ii,vv+1)     = sum(data_ym_zm_low*dlat_low)/sum(dlat_low) - sum(data_ym_zm_high*dlat_high)/sum(dlat_high) ; weighted MTG
  ;                means_deepmip(nn,ii,vv+1)     = avg(data_ym_zm_low) - avg(data_ym_zm_high) ; unweighted MTG

                  delete([/data_ym_zm_low,data_ym_zm_high,dlat_low,dlat_high/])         
                end if
            delete([/data,data_ym_zm/])         
          end if
        end do
      end do
    end do
    means_deepmip(:,:,0)  = means_deepmip(:,:,0) - 273.15
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
    means_ppe     = new((/dimsizes(models_ppe),dimsizes(exps_ppe),dimsizes(vars_ppe)+1/),"float")
    do nn=0,dimsizes(models_ppe)-1
      print("processing "+models_ppe(nn))
      do ii=0,dimsizes(exps_ppe)-1
        do vv=0,dimsizes(vars_ppe)-1
          ;--- check whether file for model+exp combination is available
          f := systemfunc("ls -r "+ppe_dir+"/"+exps_ppe(ii)+models_ppe(nn)+"*"+varnames_ppe(vv)+"*.mean.zm.nc* 2> /dev/null")
          if (.not. all(ismissing(f))) then ; file(s) found
            fi            = addfile(f,"r")
            data_ym_zm    = fi->$vars_ppe(vv)$(0,0,:,0)
            data_ym_zm!0  = "lat"
            lat          := tofloat(data_ym_zm&lat)
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
                means_ppe(nn,ii,vv)     = sum(data_ym_zm*dlat)/sum(dlat)
              ;--- SST gradient
                if (vars_ppe(vv) .eq. "temp_ym_dpth") then
                  data_ym_zm_low  = where(abs(lat) .lt. 30., data_ym_zm,data_ym_zm@_FillValue)
                  data_ym_zm_high = where(abs(lat) .gt. 60., data_ym_zm,data_ym_zm@_FillValue)
                  dlat_low        = where(abs(lat) .lt. 30., dlat,dlat@_FillValue)
                  dlat_high       = where(abs(lat) .gt. 60., dlat,dlat@_FillValue)
                
                  means_ppe(nn,ii,vv+1)     = sum(data_ym_zm_low*dlat_low)/sum(dlat_low) - sum(data_ym_zm_high*dlat_high)/sum(dlat_high) ; weighted MTG
   ;               means_ppe(nn,ii,vv+1)     = avg(data_ym_zm_low) - avg(data_ym_zm_high) ; unweighted MTG

                  delete([/data_ym_zm_low,data_ym_zm_high,dlat_low,dlat_high/])         
                end if
            delete([/data_ym_zm/])         
          end if
        end do
      end do
    end do
    means_ppe(:,:,0)  = means_ppe(:,:,0) - 273.15

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

  pl_colors             = (/"gray","black","purple","blue","cyan","green","orange"/)
  pl_styles_deepmip     = (/4,16,5,6,9,7/)
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
  
  co2_levels    = (/1,1,2,3,4,6,9/)
  
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
        poly2(count) = gsn_add_polymarker(wks,plot(2),log10(co2_levels(ll)+0.05*co2_levels(ll)),means_ppe(nn,ll,0),plres)
        if (models_ppe(nn) .ne. "f") then
          poly2_out(count) = gsn_add_polymarker(wks,plot(2),log10(co2_levels(ll)+0.05*co2_levels(ll)),means_ppe(nn,ll,0),ploutres)
        end if
        count       = count+1
      end do
    end do
    
    proxy_area_a2  = gsn_add_polygon(wks,plot(2),log10((/3.4,9,9,3.4,3.4/)),(/26,26,32,32,26/),pgres)

    ;--- SST vs. SST MTG
    res_blank@trXMinF         = 15
    res_blank@trXMaxF         = 40
    res_blank@trYMinF         = 17
    res_blank@trYMaxF         = 28
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
        poly(count) = gsn_add_polymarker(wks,plot(1),means_deepmip(nn,ll,1),means_deepmip(nn,ll,2),plres)
        count       = count+1
        plres@gsMarkerIndex   = pl_styles_ppe_f
        poly(count) = gsn_add_polymarker(wks,plot(1),means_ppe(5,ll,1),means_ppe(5,ll,2),plres)
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
        poly2(count) = gsn_add_polymarker(wks,plot(3),means_deepmip(nn,ll,1),means_deepmip(nn,ll,2),plres)
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
        poly2(count) = gsn_add_polymarker(wks,plot(3),means_ppe(nn,ll,1),means_ppe(nn,ll,2),plres)
        if (models_ppe(nn) .ne. "f") then
          poly2_out(count) = gsn_add_polymarker(wks,plot(3),means_ppe(nn,ll,1),means_ppe(nn,ll,2),ploutres)
        end if
        count       = count+1
      end do
    end do

    proxy_area_b2  = gsn_add_polygon(wks,plot(3),(/26,32,32,26,26/),(/20,20,22,22,20/),pgres)

   names_legend = array_append_record(models_deepmip,"HadCM3B_M2.1aN_PPE_f",0)
   pl_styles_legend = array_append_record(pl_styles_deepmip,pl_styles_ppe_f,0)
   
   legend = create "Legend" legendClass wks 
   "vpXF" : 0.291
   "vpYF" : .69
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
   "lgLabelFontHeightF"  : 0.007
   "lgLabelPosition" : "Right"
   "lgLabelOffsetF" : 0.1
   "lgMonoDashIndex" : False
   "lgMarkerIndexes" : pl_styles_legend
   "lgMonoMarkerColor" : True ; Only one dash pattern - solid
   "lgMarkerColor" : "black"
   "lgMarkerSizeF" : 0.01
   "lgMarkerThicknessF" : 2
   end create 

   draw(legend)

  gsn_panel(wks,plot,(/2,2/),resP)

  wks = gsn_open_wks("pdf",plot_dir+"/PI_scatter_ppe")

  plot          := new(1,graphic)
  poly          := new(2000,graphic)
  poly2         := new(2000,graphic)
  poly2_out     := new(2000,graphic)

  count         = 0
  
  res_blank@trXMinF         = 16
  res_blank@trXMaxF         = 21
  res_blank@trYMinF         = 12
  res_blank@trYMaxF         = 17
  res_blank@tmXBMode        = "Automatic"
  res_blank@tiXAxisString  = "global mean SST"
  res_blank@tiYAxisString  = "global mean surface temperature"
  res_blank@gsnLeftString  = "DeepMIP pre-industrial simulations"
  plot(0)       = gsn_csm_blank_plot(wks,res_blank)


    plres@gsMarkerOpacityF = 1.0
    do nn=0,dimsizes(models_deepmip)-1
      plres@gsMarkerSizeF   = 0.018
      plres@gsMarkerThicknessF = 3.
      plres@gsMarkerIndex = pl_styles_deepmip(nn)
      do ll=0,0
        plres@gsMarkerColor = pl_colors(ll)
        poly2(count) = gsn_add_polymarker(wks,plot(0),means_deepmip(nn,ll,1),means_deepmip(nn,ll,0),plres)
        count       = count+1
      end do
    end do   
    plres@gsMarkerOpacityF = 1.
    ploutres@gsMarkerSizeF      = 0.013
    
    do nn=0,dimsizes(models_ppe)-1
      if (models_ppe(nn) .eq. "f") then
        plres@gsMarkerIndex = pl_styles_ppe_f
        plres@gsMarkerSizeF   = 0.018
        plres@gsMarkerThicknessF = 2.
      else
        plres@gsMarkerSizeF   = 0.012
        plres@gsMarkerThicknessF = 2.
        plres@gsMarkerIndex = pl_styles_ppe
      end if
      do ll=0,0
        plres@gsMarkerColor = pl_colors(ll)
        poly2(count) = gsn_add_polymarker(wks,plot(0),means_ppe(nn,ll,1),means_ppe(nn,ll,0),plres)
        if (models_ppe(nn) .ne. "f") then
          poly2_out(count) = gsn_add_polymarker(wks,plot(0),means_ppe(nn,ll,1),means_ppe(nn,ll,0),ploutres)
        end if
        count       = count+1
      end do
    end do

   legend2 = create "Legend" legendClass wks 
   "vpXF" : 0.63
   "vpYF" : .47
   "lgAutoManage" : False
   "vpWidthF" : 0.12
   "vpHeightF" : 0.18
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
   "lgLabelFontHeightF"  : 0.012
   "lgLabelPosition" : "Right"
   "lgLabelOffsetF" : 0.1
   "lgMonoDashIndex" : False
   "lgMarkerIndexes" : pl_styles_legend
   "lgMonoMarkerColor" : True ; Only one dash pattern - solid
   "lgMarkerColor" : "black"
   "lgMarkerSizeF" : 0.016
   "lgMarkerThicknessF" : 1.5
   end create 

   draw(legend2)
    
    gsn_panel(wks,plot,(/1,1/),resP)
