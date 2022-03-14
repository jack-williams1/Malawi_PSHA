%% OPTION TO ADD NEW SOURCE ZONE

%Add two new source zones for regions not covered by Poggi SZ
%Only needs to run once

load SZ_final
clear tmp_SZ_final

mydir  = pwd; idcs   = strfind(mydir,'/');
addpath(mydir(1:idcs(end)-1));
addpath([mydir(1:idcs(end)-1) '/gis_files']);
[S,A]=shaperead('GEM_SourceZones_WGS.shp');

tmp1_geom = [S(14).X; S(14).Y]';%Coordinates from shape file, update as necessary
tmp1_SZ_final.ZONE = horzcat(tmp1_geom(~isnan(tmp1_geom(:,1))),tmp1_geom(~isnan(tmp1_geom(:,2)),2));%Coordinates from shape file, update as necessary
tmp1_SZ_final.Group = 7;%Max Group from Poggi was 6
tmp1_SZ_final.Source = length(SZ_final)+1;
tmp1_SZ_final.a = 1.74; %see spreadsheet
tmp1_SZ_final.b = 0.8; %see spreadsheet
tmp1_SZ_final.Mwmax =7.0 ; %See Fenton et al (2006)

tmp2_geom = [S(15).X; S(15).Y]';%Coordinates from shape file, update as necessary
tmp2_SZ_final.ZONE = horzcat(tmp2_geom(~isnan(tmp2_geom(:,1))),tmp2_geom(~isnan(tmp2_geom(:,2)),2));%Coordinates from shape file, update as necessary
tmp2_SZ_final.Group = 7;%Max Group from Poggi was 6
tmp2_SZ_final.Source = length(SZ_final)+2;
tmp2_SZ_final.a = 1.93; %tmp, update if necessary
tmp2_SZ_final.b = 0.8; %tmp, update if necessary
tmp2_SZ_final.Mwmax =7.0 ; %See Fenton et al (2006)


new_SZ =  [SZ_final,tmp1_SZ_final,tmp2_SZ_final];

save('new_SZ','new_SZ');