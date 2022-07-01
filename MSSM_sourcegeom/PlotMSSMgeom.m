%% Plot MSSM polygons 

%Plots each MSSM fault source polygon
%Each coordinate within polygon represents a different fault or fault section

%Load the MSSM geometry, columns are:
%1) Lat
%2) Lon
%3) Depth
%4) Source ID

source_geom=table2array(readtable('MSSM_source_geometry.csv'));
MSSM_id=[min(source_geom(:,4)):1:max(source_geom(:,4))];

s_th = 35; % Thickness of seismogenic crust
Region1 = [32.6 36 -17.2 -9]; %Bounds of figure plot

figure(1);
for ii=1:length(MSSM_id)
    tmp_indx=find(MSSM_id(ii)==source_geom(:,4));
    
    plot3(source_geom(tmp_indx,2),source_geom(tmp_indx,1),source_geom(tmp_indx,3),'k-');hold on
   
end

xlim([Region1(1) Region1(2)]); ylim([Region1(3) Region1(4)]); grid on;
ax=gca; daspect([0.01 0.01 1]); zlim([0 s_th]); set(ax,'ZDir','Reverse'); hold on 
view(30,45)

