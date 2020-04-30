deepmip_dir = "/export/acrc/DeepMIP_Model_Output_shared"
ppe_dir     = "/home/bridge/wb19586/spinup_validation/data/fields"
plot_dir    = "/home/bridge/wb19586/analysis/ppe_deepmip/plots"

models_deepmip    = (/"HadCM3BL_M2.1aN","HadCM3B_M2.1aN","GFDL_CM2.1","CESM1.2_CAM5","COSMOS-landveg_r2413"/)
versions_deepmip  = (/"v1.0","v1.0","v1.0","v1.0","v1.0"/)
exps_deepmip      = (/"{pi_control,piControl}","1x","2x","3x","4x","6x","9x"/)
vars_deepmip      = (/"tas","tos"/)

;--- load DeepMIP TAS+TOS and calculate means
  means_deepmip     = new((/dimsizes(models_deepmip),dimsizes(exps_deepmip),dimsizes(vars_deepmip)+1/),"float")
  do nn=0,dimsizes(models_deepmip)-1
    do ii=0,dimsizes(exps_deepmip)-1
      do vv=0,dimsizes(vars_deepmip)-1
        ;--- check whether file for model+exp combination is available
        f := systemfunc("ls "+deepmip_dir+"/*/"+models_deepmip(nn)+"/*"+exps_deepmip(ii)+"*/"+versions_deepmip(nn)+"/*"+vars_deepmip(vv)+"-"+versions_deepmip(nn)+".mean.* 2> /dev/null")
        if (.not. all(ismissing(f))) then ; file(s) found
          if (dimsizes(f) .eq. 1) then    ; 1 file found
            fi        = addfile(f,"r")
            data      = fi->$vars_deepmip(vv)$
          else if (dimsizes(f) .eq. 2) then    ; original+remapped grid for "tos" data available
            fi         = addfile(f(0),"r")
            data       = fi->$vars_deepmip(vv)$
            data_dims := dimsizes(rm_single_dims(data))
            if (.not. any(data_dims .eq. 360)) then ; not 1x1 degree grid, therefore load second file
              delete(data)
              fi         = addfile(f(1),"r")
              data       = fi->$vars_deepmip(vv)$
              data_dims := dimsizes(rm_single_dims(data))
              if (.not. any(data_dims .eq. 360)) then  ; still not 1x1 degree grid -> unexpected
                print("no regular 1x1 degree grid found. EXIT!")
                print("error for "+models_deepmip(nn)+"/"+exps_deepmip(ii)+"/"+vars_deepmip(vv))
                exit
              end if
            end if
          end if
          end if
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

                means_deepmip(nn,ii,vv+1)     = sum(data_ym_zm_low*dlat_low)/sum(dlat_low) - sum(data_ym_zm_high*dlat_high)/sum(dlat_high)
                delete([/data_ym_zm_low,data_ym_zm_high,dlat_low,dlat_high/])         
              end if
          delete([/data,data_ym_zm/])         
        end if
      end do
    end do
  end do

means_deepmip(:,:,0)  = means_deepmip(:,:,0) - 273.15

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
  plres                 = True
  plres@gsMarkerSizeF   = 0.014
  plres@gsMarkerThicknessF = 2.
  pl_colors             = (/"gray","black","purple","blue","cyan","green","orange"/)
  pl_styles             = (/4,16,5,6,9/)
  
  pgres                 = True
  pgres@gsFillColor     = "transparent"      ;-- fill color
  pgres@gsEdgesOn       = True
  pgres@gsEdgeColor     = "red"         ;-- set opacity of polygon
  pgres@gsEdgeThicknessF = 2.0


  plot          = new(2,graphic)
  poly          = new(200,graphic)
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
      plres@gsMarkerIndex = pl_styles(nn)
      do ll=0,dimsizes(co2_levels)-1
        plres@gsMarkerColor = pl_colors(ll)
        poly(count) = gsn_add_polymarker(wks,plot(0),log10(co2_levels(ll)),means_deepmip(nn,ll,0),plres)
        count       = count+1
      end do
    end do
    proxy_area_a  = gsn_add_polygon(wks,plot(0),log10((/3.3,9,9,3.3,3.3/)),(/26,26,32,32,26/),pgres)

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
      plres@gsMarkerIndex = pl_styles(nn)
      do ll=0,dimsizes(co2_levels)-1
        plres@gsMarkerColor = pl_colors(ll)
        poly(count) = gsn_add_polymarker(wks,plot(1),means_deepmip(nn,ll,1),means_deepmip(nn,ll,2),plres)
        count       = count+1
      end do
    end do
    proxy_area_b  = gsn_add_polygon(wks,plot(1),(/26,32,32,26,26/),(/20,20,22,22,20/),pgres)

   names = models_deepmip

   legend = create "Legend" legendClass wks 
   "vpXF" : 0.59
   "vpYF" : .756
   "lgAutoManage" : False
   "vpWidthF" : 0.08
   "vpHeightF" : 0.15
   "lgLineThicknessF" : 8.0
   "lgLineDashSegLenF" : 0.08
   "lgPerimOn" : True
   "lgLabelJust" : "CenterLeft"
   "lgItemCount" : dimsizes(names)
   "lgItemOrder" : ispan(dimsizes(names)-1,0,1)
   "lgLabelStrings" : names
   "lgItemType" : "Markers"
   "lgLabelsOn" : True 
   "lgLineLabelsOn" : False 
   "lgLabelFontHeightF"  : 0.01
   "lgLabelPosition" : "Right"
   "lgLabelOffsetF" : 0.1
   "lgMonoDashIndex" : False
   "lgMarkerIndexes" : pl_styles
   "lgMonoMarkerColor" : True ; Only one dash pattern - solid
   "lgMarkerColor" : "black"
   "lgMarkerSizeF" : 0.014
   "lgMarkerThicknessF" : 2
   end create 

   draw(legend)

  gsn_panel(wks,plot,(/2,1/),resP)
    