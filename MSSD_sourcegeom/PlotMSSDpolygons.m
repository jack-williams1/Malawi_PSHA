%% Plot MSSD polygons with the MAFD

%Plots each MSSD fault source polygon
%Each coordinate within polygon represents a different fault section

%Note requires Matlab mapping toolbox to convert MAFD coordinates to wgs84

%Load the MSSD geometry, columns are:
%1) Lat
%2) Lon
%3) Depth
%4) Source ID

mydir  = pwd; idcs   = strfind(mydir,'/');
addpath(mydir(1:idcs(end)-1));
addpath([mydir(1:idcs(end)-1) '/gis_files']);

load fault_geom_1; num_fault=length(cellfun(@numel,fault_geometry_points_1));

%Read in MSSD shape file, update if necessary
MSSD = shaperead('MSSD_fault.shp');
num_traces = length(MSSD);

ellipsoid = almanac('earth','wgs84','meters');
utms = defaultm('utm'); utms.zone = '36L'; utms.geoid = ellipsoid; 
utms.flatlimit = []; utms.maplatlimit = []; utms = defaultm(utms);

s_th = 35; % Thickness of seismogenic crust
Region1 = [32.6 36 -17.2 -9]; %Bounds of figure plot

figure(3);

%Plot MSSD polygons
for i=1:num_fault
  plot3(fault_geometry_points_1{i}(:,2),fault_geometry_points_1{i}(:,1),fault_geometry_points_1{i}(:,3),'k-'); 
  hold on
end

hold on

%Plot MSSD traces
for kk = 1:num_traces    
    [LAT,LON] = minvtran(utms,MSSD(kk).X,MSSD(kk).Y);
    plot3(LON,LAT,ones(length(LON),1),'r-','linewidth',2); hold on; 
end

xlim([Region1(1) Region1(2)]); ylim([Region1(3) Region1(4)]); grid on;
ax=gca; daspect([0.01 0.01 1]); zlim([0 s_th]); set(ax,'ZDir','Reverse'); hold on 
view(30,45)

