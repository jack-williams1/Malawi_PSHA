%% Plot Hodge 2015 PSHA map for Malawi for 2%  PoE in 50 years
% Then also plots comparison with SAFER-PREPARED PSHA map

mydir  = pwd; idcs   = strfind(mydir,'/');
addpath(mydir(1:idcs(end)-1));
addpath(mydir(1:idcs(end-1)-1));
addpath([mydir(1:idcs(end-1)-1) '/misc_functions']);
addpath([mydir(1:idcs(end-1)-1) '/gis_files']);

load('syncat_PSHA_MSSM_input','Region','site_grid_interval','vs30_site_ref','prob_level')
prob_level=prob_level(2); %Hodge et al 2015 only compared for 2% PoE in 50 years

% Load PSHA files from Hodge 2015
load h2015_map_values
  

%% Create new grid for MSSM as grid interval is different from Hodge et al 2015

%find overlapping region between MSSM and Hodge2015 maps
lat_range = [max(min(site(:,1)),Region(3)) min(max(site(:,1)),Region(4))];
lon_range = [max(min(site(:,2)),Region(1)) min(max(site(:,2)),Region(2))];    

MSSM_malawi_lat_site_vec = lat_range(1):site_grid_interval:lat_range(2); %can be changed for more accurate plotting
MSSM_malawi_lon_site_vec = lon_range(1):site_grid_interval:lon_range(2);
[Lat_site,Lon_site] = meshgrid(MSSM_malawi_lat_site_vec,MSSM_malawi_lon_site_vec);
MSSM_site = [Lat_site(:) Lon_site(:)];

for pp=1:length(prob_level)

clear pga_mesh
    for ii=1:length(MSSM_site)
         [~,mindiff_id] = min((MSSM_site(ii,1)-site(:,1)).^2 + (MSSM_site(ii,2)-site(:,2)).^2);
       %Find closes PSHA values to those in MSSM map
         pga(ii,:) = PSA_fractile_CDF(:,mindiff_id);
         
    end

pga_mesh = reshape(pga,[length(MSSM_malawi_lon_site_vec),length(MSSM_malawi_lat_site_vec)]);
pga_meshh{pp} = rot90(pga_mesh,3);
pga_meshh{pp} = flip(pga_meshh{pp},2);

end

%% Load map data

ellipsoid = almanac('earth','wgs84','meters');
utms = defaultm('utm'); utms.zone = '36L'; utms.geoid = ellipsoid; 
utms.flatlimit = []; utms.maplatlimit = []; utms = defaultm(utms);

%Read in seismogenic sources geographic extent
MSSM = shaperead('MSSM_Faults.shp');
h_sources = shaperead('hodge2015faultsources.shp');
h_sources(1).Name='BMF';
num_traces = length(MSSM);

load map_data_EastAfrica.mat;
LakeMalawi = shaperead('malawi_lake.shp');
LakeMalawiCoord = [LakeMalawi.Y(1,1:end-1)' LakeMalawi.X(1,1:end-1)'];

%Adapt figure axis to just around lake
Region_2 =  [lon_range(1) lon_range(2) -15 lat_range(2)];
h2015_vs30 = 760; 

plabel_opt=strcat(["2% PoE in 50 years"]);
label_opt=strcat(["(a) H2015","(b) MSSM","(c) MSSM-H2015"]);


%% Compare to SAFER-PREPARED PSHA map (USGS vs30 and vs30 ref)

for pp=1:length(prob_level)

s_prepare_pga=readtable(strcat('MalawiPGA_USGS_vs30_',string(prob_level(pp)))); s_prepare_pga_ref=readtable(strcat('MalawiPGA_ref_vs30_',string(prob_level(pp))));
s_prepare_pga=table2array(s_prepare_pga); s_prepare_pga_ref=table2array(s_prepare_pga_ref);
s_prepare_pga_array=s_prepare_pga(find(s_prepare_pga(:,1) >= lat_range(1) & s_prepare_pga(:,1) <= lat_range(2)),find(s_prepare_pga(1,:) >= lon_range(1) & s_prepare_pga(1,:) <= lon_range(2)));
s_prepare_pga_array_ref=s_prepare_pga_ref(find(s_prepare_pga_ref(:,1) >= lat_range(1) & s_prepare_pga_ref(:,1) <= lat_range(2)),find(s_prepare_pga_ref(1,:) >= lon_range(1) & s_prepare_pga_ref(1,:) <= lon_range(2)));

s_prepare_pga_array1{pp} = flipud(s_prepare_pga_array);
s_prepare_pga_array1_ref{pp} = flipud(s_prepare_pga_array_ref);

pga_ratio{pp}=zeros(height(s_prepare_pga_array),width(s_prepare_pga_array));
pga_ratio_ref{pp}=zeros(height(s_prepare_pga_array_ref),width(s_prepare_pga_array_ref));

for i=1:height(s_prepare_pga_array)
    for j=1:width(s_prepare_pga_array)
        pga_ratio{pp}(i,j)=s_prepare_pga_array1{pp}(i,j)-pga_meshh{pp}(i,j);
        pga_ratio_ref{pp}(i,j)=s_prepare_pga_array1_ref{pp}(i,j)-pga_meshh{pp}(i,j);
    end
end

end
     
%% Plot figure

cmap = crameri('batlow');

vs30_opt=2; %set to 1 for USGS vs30, 2 for ref vs_30

if vs30_opt==1
    
  for i=1:length(prob_level)  
    tmp1{pp}=s_prepare_pga_array1{pp};
    tmp2{pp}=pga_ratio{pp};
    vs_val='USGS vs30';
  end 
  
else    
    
   for i=1:length(prob_level)  
    tmp1{pp}=s_prepare_pga_array1_ref{pp};
    tmp2{pp}=pga_ratio_ref{pp};
    vs_val='vs30 = 760 m/s';
   end 
  
end

figure(2);
tiledlayout(length(prob_level),3,'tilespacing','compact')


for pp=1:length(prob_level)

nexttile
surf(YI(:,1)',XI(1,:)',flip(rot90(ZI_2475yrs,3),2)); hold on; shading interp; view(0,90); axis equal; axis(Region_2); title(label_opt(pp,1)); subtitle({plabel_opt(pp),['vs30 = ' num2str(h2015_vs30) ' m/s']});  set(gca,'fontsize',13.5);        
plot3(MapData2(:,2),MapData2(:,1),1000*ones(length(MapData2(:,1)),1),'w-');
fill3(LakeMalawiCoord(:,2),LakeMalawiCoord(:,1),1000*ones(length(LakeMalawiCoord),1),'w');
hold on; colormap(gca,cmap); caxis([0.0 0.8]); h1=colorbar; h1.Label.String = 'PGA (g)';

        t_corr_x=[-0.4 0 0 -1.1 0.1 0 0.4];
        t_corr_y=[-0.5 0 0 0 0 0 -0.3];

        for jj = 1:length(h_sources)    
        plot3(h_sources(jj).X,h_sources(jj).Y,1000*ones(length(h_sources(jj).Y),1),'r','linewidth',2); hold on;

            if jj==1 || jj==4 || jj==5 || jj ==7 
                text(h_sources(jj).X(1)+t_corr_x(jj),h_sources(jj).Y(1)+t_corr_y(jj),1000,h_sources(jj).Name,'Color', 'white','FontSize',10.5);
            end
        end

hold off

nexttile
surf(MSSM_malawi_lon_site_vec,MSSM_malawi_lat_site_vec,tmp1{pp}); hold on; shading interp; view(0,90); axis equal; axis(Region_2); title(label_opt(pp,2)); subtitle([plabel_opt(pp),vs_val]); set(gca,'fontsize',13.5);         
plot3(MapData2(:,2),MapData2(:,1),1000*ones(length(MapData2(:,1)),1),'w-');
fill3(LakeMalawiCoord(:,2),LakeMalawiCoord(:,1),1000*ones(length(LakeMalawiCoord),1),'w');
hold on;  colormap(gca,cmap); caxis([0.0 0.8]); h2=colorbar; h2.Label.String = 'PGA (g)';
        for jj = 1:num_traces    
            plot3(MSSM(jj).X,MSSM(jj).Y,1000*ones(length(MSSM(jj).X),1),'r','linewidth',0.75); hold on;  
        end
      
hold off        
nexttile
cmap1 = crameri('vik'); 
surf(MSSM_malawi_lon_site_vec,MSSM_malawi_lat_site_vec,tmp2{pp}); hold on; shading interp; view(0,90); axis equal; axis(Region_2); title([label_opt(pp,3),'']);  set(gca,'fontsize',13.5);        
plot3(MapData2(:,2),MapData2(:,1),1000*ones(length(MapData2(:,1)),1),'k-');
fill3(LakeMalawiCoord(:,2),LakeMalawiCoord(:,1),1000*ones(length(LakeMalawiCoord),1),'w');
hold on;  colormap(gca,cmap1);  h3=colorbar; h3.Label.String = 'PGA difference (g)'; caxis([-0.4 0.4]);
      %{
      for jj = 1:num_traces    
        plot3(MSSM(jj).X,MSSM(jj).Y,1000*ones(length(MSSM(jj).X),1),'r','linewidth',0.75); hold on;  
      end
      %}
cmap = crameri('batlow');
  
end    
set(gcf,'position',[785 84 838 455]) 

%% Save data for combined plot with GEM comparison in gem_psha_comparison  

h2015_pga_ratio_ref=pga_ratio_ref; h2015_pga_ratio=pga_ratio;
h15_pga_meshh=pga_meshh; h15_prob_level=prob_level;
h15_MSSM_malawi_lat_site_vec=MSSM_malawi_lat_site_vec;
h15_MSSM_malawi_lon_site_vec=MSSM_malawi_lon_site_vec;

save('h2015_map_pga','XI','YI','ZI_2475yrs','h2015_pga_ratio','h2015_pga_ratio_ref',...
    'h15_prob_level','Region_2','h2015_vs30','h_sources','t_corr_x','t_corr_y',...
    'h15_MSSM_malawi_lat_site_vec','h15_MSSM_malawi_lon_site_vec');