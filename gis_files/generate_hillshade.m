%% Generate hillshade from DEM for Malawi

%Use in source geomfigure
%Requires Matlab Image Processing Toolbox

load Malawi_GDEM2 %DEM

mydir  = pwd; idcs   = strfind(mydir,'/');
addpath(mydir(1:idcs(end)-1));
addpath([mydir(1:idcs(end)-1) '/misc_functions']);

inputexcel = 'syncat_PSHA_MSSD_input.xlsx';

ellipsoid = almanac('earth','wgs84','meters');
utms = defaultm('utm'); utms.zone = '36L'; utms.geoid = ellipsoid; 
utms.flatlimit = []; utms.maplatlimit = []; utms = defaultm(utms);

%Region bound by hillshade
Region1 = [xlsread(inputexcel,1,'D47') xlsread(inputexcel,1,'F47') xlsread(inputexcel,1,'C47') xlsread(inputexcel,1,'E47')];

%% Clean and downsample DEM

clear tmp
extra = 50;
[~,tmp(1)] = min(abs(Lat_Malawi-Region1(4))); tmp(1) = max(tmp(1)-extra,1);
[~,tmp(2)] = min(abs(Lat_Malawi-Region1(3))); tmp(2) = tmp(2) + extra;
[~,tmp(3)] = min(abs(Lon_Malawi-Region1(1))); tmp(3) = max(tmp(3)-extra,1);
[~,tmp(4)] = min(abs(Lon_Malawi-Region1(2))); tmp(4) = tmp(4);

resamp_interval = 15;
Lat_Malawi   = Lat_Malawi(tmp(1):resamp_interval:tmp(2));
Lon_Malawi   = Lon_Malawi(tmp(3):resamp_interval:tmp(4));
GDEM2_Malawi = GDEM2_Malawi(tmp(1):resamp_interval:tmp(2),tmp(3):resamp_interval:tmp(4));
save('hillshade_xy_data','Lon_Malawi','Lat_Malawi')
%% Generate hillshade

%z factor set for difference between metre and lat and lon units
%appropriate for 10 degree lat
h=hillshade(GDEM2_Malawi,Lon_Malawi,Lat_Malawi,'zfactor',0.00000912); 
imwrite(ind2rgb(im2uint8(mat2gray(h)), gray(256)), 'malawi_hillshade.jpg');
