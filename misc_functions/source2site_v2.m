% RANDOMLY FLOAT RUPTURES ON FAULT PLANE

% Associate ruptures with random grid point from faultgrid_wgs84
% Ruptures from Adapted MSSD event catalog 
% For events with Mw<5.4

function [gp]= source2site_v2(syncat_ijk,faultgrid_wgs84_1,faultgrid_wgs84_2,grid_indx1,grid_indx2)

  clear gp
  
    if syncat_ijk(5)==1  %Partial fault width, use fault_geom_1 model
        grd_indx = find(syncat_ijk(4)==grid_indx1(:));
        r_indx = randsample([1:1:length(faultgrid_wgs84_1{grd_indx})],1);
        gp=faultgrid_wgs84_1{grd_indx}(r_indx,1:3);

    else %Full crust fault width, use fault_geom_2 model
        grd_indx = find(syncat_ijk(4)==grid_indx2(:));
        r_indx = randsample([1:1:length(faultgrid_wgs84_2{grd_indx})],1);
        gp=faultgrid_wgs84_2{grd_indx}(r_indx,1:3);
    end
    
end
