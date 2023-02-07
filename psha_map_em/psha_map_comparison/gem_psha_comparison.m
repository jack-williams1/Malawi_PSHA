%% Plots GEM PSHA map for Malawi for 10% PoE

% Then also plots comparison with SAFER-PREPARED PSHA map
% and loads in Hodge 2015 to compare that map too
mydir  = pwd; idcs   = strfind(mydir,'/');
addpath(mydir(1:idcs(end)-1)); addpath(mydir(1:idcs(end-1)-1));
addpath([mydir(1:idcs(end-1)-1) '/misc_functions']);
addpath([mydir(1:idcs(end-1)-1) '/gis_files']);

load('syncat_PSHA_MSSM_input','Region','site_grid_interval','vs30_site_ref','prob_level')

malawi_Lon_w=Region(1); malawi_Lon_e=Region(2);
malawi_Lat_s=Region(3); malawi_Lat_n=Region(4);

load gem_psha_Malawi %load GEM data

gem_vs30 = 600; %See Poggi et al 2017

%load H2015 pga comparisons
load h2015_map_pga; 
    
gem_psha_Malawi_array{1}=table2array(gem_psha_Malawi(:,[1:3])); % 10% PoE in 50 years
gem_psha_Malawi_array{2}=table2array(gem_psha_Malawi(:,[1:2,4])); % 2% PoE in 50 years

%% Create grid of values

malawi_lat_site_vec = malawi_Lat_s:site_grid_interval:malawi_Lat_n; %can be changed for more accurate plotting
malawi_lon_site_vec = malawi_Lon_w:site_grid_interval:malawi_Lon_e;
[Lat_site,Lon_site] = meshgrid(malawi_lat_site_vec,malawi_lon_site_vec);
site = [Lat_site(:) Lon_site(:)];

pga_meshh={};

for pp=1:length(prob_level)

    for ii = 1:length(site)
        [~,mindiff_id] = min((gem_psha_Malawi_array{pp}(:,2)-site(ii,1)).^2 + (gem_psha_Malawi_array{pp}(:,1)-site(ii,2)).^2);
        %Find closes PSHA values to those in PREPARE map
        pga(ii,1) = gem_psha_Malawi_array{pp}(mindiff_id,3);
    end


pga_mesh = reshape(pga,[length(malawi_lon_site_vec),length(malawi_lat_site_vec)]);
pga_meshh{pp} = rot90(pga_mesh,3);
pga_meshh{pp} = flip(pga_meshh{pp},2);
end

%% Load map data

load map_data_EastAfrica.mat;
ellipsoid = almanac('earth','wgs84','meters');
utms = defaultm('utm'); utms.zone = '36L'; utms.geoid = ellipsoid; 
utms.flatlimit = []; utms.maplatlimit = []; utms = defaultm(utms);
%Read in MSSM shape file, update if necessary
MSSM = shaperead('MSSM_Faults.shp');
num_traces = length(MSSM);

LakeMalawi = shaperead('malawi_lake.shp');
LakeMalawiCoord = [LakeMalawi.Y(1,1:end-1)' LakeMalawi.X(1,1:end-1)'];
load LakeMalawiBorder

label_opt=vertcat(strcat(["(a) P2017","(b) MSSM","(c) MSSM-P2017"]),strcat(["(d) H2015","(e) MSSM","(f) MSSM-H2015"]));



%% Compare to SAFER-PREPARED PSHA map (USGS vs30 and vs30 ref)

for pp=1:length(prob_level)

s_prepare_pga=readtable(strcat('MalawiPGA_USGS_vs30_',string(prob_level(pp)))); s_prepare_pga_ref=readtable(strcat('MalawiPGA_ref_vs30_',string(prob_level(pp))));

s_prepare_pga=s_prepare_pga(2:height(s_prepare_pga),2:width(s_prepare_pga));%remove lat and lon column headers
s_prepare_pga_ref=s_prepare_pga_ref(2:height(s_prepare_pga_ref),2:width(s_prepare_pga_ref));%remove lat and lon column headers

s_prepare_pga_array=table2array(s_prepare_pga);
s_prepare_pga_array_ref=table2array(s_prepare_pga_ref);

s_prepare_pga_array1{pp} = flipud(s_prepare_pga_array);
s_prepare_pga_array1_ref{pp} = flipud(s_prepare_pga_array_ref);

pga_ratio{pp}=zeros(height(s_prepare_pga_array),width(s_prepare_pga));
pga_ratio_ref{pp}=zeros(height(s_prepare_pga_array_ref),width(s_prepare_pga_ref));

for i=1:height(s_prepare_pga)
    for j=1:width(s_prepare_pga)
        pga_ratio{pp}(i,j)=s_prepare_pga_array1{pp}(i,j)-pga_meshh{pp}(i,j);
        pga_ratio_ref{pp}(i,j)=s_prepare_pga_array1_ref{pp}(i,j)-pga_meshh{pp}(i,j);
    end
end

end
      
%% Plot figure

vs30_opt=2; %set to 1 for USGS vs30, 2 for ref vs_30

if vs30_opt==1
    
  for pp=1:length(prob_level)  
    tmp1{pp}=s_prepare_pga_array1{pp};
    tmp2{pp}=pga_ratio{pp}; 
    vs_val='USGS vs30';
  end 
  
  h15_tmp2=h2015_pga_ratio;
  
else    
    
   for pp=1:length(prob_level)  
    tmp1{pp}=s_prepare_pga_array1_ref{pp};
    tmp2{pp}=pga_ratio_ref{pp};
    vs_val='vs30 = 760 m/s';
   end 
  
   h15_tmp2=h2015_pga_ratio_ref;
  
end


%% Plot combined GEM (10% PoE in 50 years) and H2015 comparison (2% PoE in 50 years)

%Select values prob levels for comparison
%1) 10% PoE in 50 years/500 year return period
%2) 2% PoE in 50 years
pp = 1;

if pp==1  
    ZI=ZI_500yrs; c_axis_u=0.4; plabel_opt=strcat(["10% PoE in 50 years","500 yr return period"]);
elseif pp==2
    ZI=ZI_2475yrs; c_axis_u=0.8; plabel_opt=strcat(["2% PoE in 50 years","2% PoE in 50 years"]);
end

cmap = crameri('batlow');

figure(1);
tiledlayout(2,3,'tilespacing','tight')

%(a) P2017
nexttile
surf(malawi_lon_site_vec,malawi_lat_site_vec,pga_meshh{pp}); hold on; shading interp; view(0,90); axis equal; axis(Region); title(label_opt(1,1),'fontweight','normal'); subtitle({plabel_opt(1),['vs30 = ' num2str(gem_vs30) ' m/s']});       

plot3(MapData2(:,2),MapData2(:,1),1000*ones(length(MapData2(:,1)),1),'w-');
%plotting outline of Lake Malawi, not its fill!
plot3(LakeMalawiCoord(:,2),LakeMalawiCoord(:,1),1000*ones(length(LakeMalawiCoord),1),'b-');
%plot3(LakeMalawiBorder(:,2),LakeMalawiBorder(:,1),1000*ones(length(LakeMalawiBorder),1),'k-');hold on
hold on; colormap(gca,cmap); caxis([0.0 c_axis_u]); set(gca,'fontsize',11); 
 
%(b) MSSM
nexttile
surf(malawi_lon_site_vec,malawi_lat_site_vec,tmp1{pp}); hold on; shading interp; view(0,90); axis equal; axis(Region); title(label_opt(1,2),'fontweight','normal'); subtitle([plabel_opt(1),vs_val]);        
plot3(MapData2(:,2),MapData2(:,1),1000*ones(length(MapData2(:,1)),1),'w-');
%plotting outline of Lake Malawi, not its fill!
plot3(LakeMalawiCoord(:,2),LakeMalawiCoord(:,1),1000*ones(length(LakeMalawiCoord),1),'b-');
%plot3(LakeMalawiBorder(:,2),LakeMalawiBorder(:,1),1000*ones(length(LakeMalawiBorder),1),'k-');hold on
hold on;  colormap(gca,cmap); caxis([0.0 c_axis_u]); h2=colorbar; h2.Label.String = 'PGA (g)';
set(gca,'XTick',[], 'YTick', [],'fontsize',11,'TitleFontSizeMultiplier',1.09);
        for jj = 1:num_traces    
        plot3(MSSM(jj).X,MSSM(jj).Y,1000*ones(length(MSSM(jj).X),1),'r','linewidth',0.75); hold on;  
        end
hold off

%(c) MSSM-P2017               
nexttile
cmap1 = crameri('vik'); 
surf(malawi_lon_site_vec,malawi_lat_site_vec,tmp2{pp}); hold on; shading interp; view(0,90); axis equal; axis(Region); title([label_opt(1,3),''],'fontweight','normal');         
plot3(MapData2(:,2),MapData2(:,1),1000*ones(length(MapData2(:,1)),1),'k-');
%plotting outline of Lake Malawi, not its fill!
plot3(LakeMalawiCoord(:,2),LakeMalawiCoord(:,1),1000*ones(length(LakeMalawiCoord),1),'b-');
hold on;  colormap(gca,cmap1);  h3=colorbar; h3.Label.String = 'PGA difference (g)'; caxis([-0.4 0.4]);
set(gca,'XTick',[], 'YTick', [],'fontsize',11,'TitleFontSizeMultiplier',1.09);
      for jj = 1:num_traces    
        plot3(MSSM(jj).X,MSSM(jj).Y,1000*ones(length(MSSM(jj).X),1),'r','linewidth',0.75); hold on;  
      end
      

cmap = crameri('batlow');

%(d) H2015
nexttile
surf(YI(:,1)',XI(1,:)',flip(rot90(ZI,3),2)); hold on; shading interp; view(0,90); axis equal; axis(Region_2); title(label_opt(2,1),'fontweight','normal'); subtitle({plabel_opt(2),['vs30 = ' num2str(h2015_vs30) ' m/s']});      
plot3(MapData2(:,2),MapData2(:,1),1000*ones(length(MapData2(:,1)),1),'w-');
%plotting outline of Lake Malawi, not its fill!
plot3(LakeMalawiCoord(:,2),LakeMalawiCoord(:,1),1000*ones(length(LakeMalawiCoord),1),'b-');
hold on; colormap(gca,cmap); caxis([0.0 c_axis_u]); set(gca,'fontsize',11); 

        t_corr_x=[-0.4 0 0 -1.1 0.1 0 0.4];
        t_corr_y=[-0.5 0 0 0 0 0 -0.3];

        for jj = 1:length(h_sources)    
        plot3(h_sources(jj).X,h_sources(jj).Y,1000*ones(length(h_sources(jj).Y),1),'r','linewidth',2); hold on;

            if jj==1 || jj==4 || jj==5 || jj ==7 
                text(h_sources(jj).X(1)+t_corr_x(jj),h_sources(jj).Y(1)+t_corr_y(jj),1000,h_sources(jj).Name,'Color', 'white','FontSize',10.5);
            end
        end

hold off

%(e) MSSM
nexttile
surf(malawi_lon_site_vec,malawi_lat_site_vec,tmp1{pp}); hold on; shading interp; view(0,90); axis equal; axis(Region_2); title(label_opt(2,2),'fontweight','normal'); subtitle([plabel_opt(2),vs_val]);      
plot3(MapData2(:,2),MapData2(:,1),1000*ones(length(MapData2(:,1)),1),'w-');
%plotting outline of Lake Malawi, not its fill!
plot3(LakeMalawiCoord(:,2),LakeMalawiCoord(:,1),1000*ones(length(LakeMalawiCoord),1),'b-');
hold on;  colormap(gca,cmap); caxis([0.0 c_axis_u]); h2=colorbar; h2.Label.String = 'PGA (g)';
set(gca,'XTick',[], 'YTick', [],'fontsize',11,'TitleFontSizeMultiplier',1.09); 
        for jj = 1:num_traces    
            plot3(MSSM(jj).X,MSSM(jj).Y,1000*ones(length(MSSM(jj).X),1),'r','linewidth',0.75); hold on;  
        end
        
hold off  

%(f) MSSM-H2015
nexttile
cmap1 = crameri('vik'); 
surf(h15_MSSM_malawi_lon_site_vec,h15_MSSM_malawi_lat_site_vec,h15_tmp2{pp}); hold on; shading interp; view(0,90); axis equal; axis(Region_2); title([label_opt(2,3),''],'fontweight','normal');      
plot3(MapData2(:,2),MapData2(:,1),1000*ones(length(MapData2(:,1)),1),'k-');
%plotting outline of Lake Malawi, not its fill!
plot3(LakeMalawiCoord(:,2),LakeMalawiCoord(:,1),1000*ones(length(LakeMalawiCoord),1),'b-');
hold on;  colormap(gca,cmap1);  h3=colorbar; h3.Label.String = 'PGA difference (g)'; caxis([-0.4 0.4]);
set(gca,'XTick',[], 'YTick', [],'fontsize',11,'TitleFontSizeMultiplier',1.09); 
      %{
      for jj = 1:num_traces    
        plot3(MSSM(jj).X,MSSM(jj).Y,1000*ones(length(MSSM(jj).X),1),'r','linewidth',0.75); hold on;  
      %}

  
set(gcf,'position',[785 84 838 895])         

%% Plot GEM comparison for both prob_levels

cmap = crameri('batlow');
figure(2);
tiledlayout(2,3,'tilespacing','tight')

for pp=1:length(prob_level)

nexttile
surf(malawi_lon_site_vec,malawi_lat_site_vec,pga_meshh{pp}); hold on; shading interp; view(0,90); axis equal; axis(Region); title(label_opt(pp,1)); subtitle({plabel_opt(pp),['vs30 = ' num2str(gem_vs30) ' m/s']});  set(gca,'fontsize',13.5);        
plot3(MapData2(:,2),MapData2(:,1),1000*ones(length(MapData2(:,1)),1),'w-');
fill3(LakeMalawiCoord(:,2),LakeMalawiCoord(:,1),1000*ones(length(LakeMalawiCoord),1),'w');
hold on; colormap(gca,cmap); caxis([0.0 1]); h1=colorbar; h1.Label.String = 'PGA (g)';
 
hold off

nexttile
surf(malawi_lon_site_vec,malawi_lat_site_vec,tmp1{pp}); hold on; shading interp; view(0,90); axis equal; axis(Region); title(label_opt(pp,2)); subtitle([plabel_opt(pp),vs_val]); set(gca,'fontsize',13.5);         
plot3(MapData2(:,2),MapData2(:,1),1000*ones(length(MapData2(:,1)),1),'w-');
fill3(LakeMalawiCoord(:,2),LakeMalawiCoord(:,1),1000*ones(length(LakeMalawiCoord),1),'w');
hold on;  colormap(gca,cmap); caxis([0.0 1]); h2=colorbar; h2.Label.String = 'PGA (g)';
        for jj = 1:num_traces    
            plot3(MSSM(jj).X,MSSM(jj).Y,1000*ones(length(MSSM(jj).X),1),'r','linewidth',0.75); hold on;  
        end
      
hold off        
nexttile
cmap1 = crameri('vik'); 
surf(malawi_lon_site_vec,malawi_lat_site_vec,tmp2{pp}); hold on; shading interp; view(0,90); axis equal; axis(Region); title([label_opt(pp,3),'']);  set(gca,'fontsize',13.5);        
plot3(MapData2(:,2),MapData2(:,1),1000*ones(length(MapData2(:,1)),1),'k-');
fill3(LakeMalawiCoord(:,2),LakeMalawiCoord(:,1),1000*ones(length(LakeMalawiCoord),1),'w');
hold on;  colormap(gca,cmap1);  h3=colorbar; h3.Label.String = 'PGA difference (g)'; caxis([-0.4 0.4]);
      for jj = 1:num_traces    
            plot3(MSSM(jj).X,MSSM(jj).Y,1000*ones(length(MSSM(jj).X),1),'r','linewidth',0.75); hold on;  
      end
      
cmap = crameri('batlow');
  
end    
set(gcf,'position',[785 84 838 895])      
      

%% PGA ratio analysis
 
 for pp=1:length(prob_level)

gem_analysis(pp,1:2)=[median(gem_psha_Malawi_array{pp}(:,3)) mean(gem_psha_Malawi_array{pp}(:,3))];    
     
MSSM_pga=reshape(s_prepare_pga_array1{pp},[length(malawi_lon_site_vec)*length(malawi_lat_site_vec),1]);  
MSSM_analysis(pp,1:3) =  [median(MSSM_pga) mean(MSSM_pga) max(MSSM_pga)];

MSSM_pga_ref=reshape(s_prepare_pga_array1_ref{pp},[length(malawi_lon_site_vec)*length(malawi_lat_site_vec),1]);  
MSSM_analysis_ref(pp,1:3) =  [median(MSSM_pga_ref) mean(MSSM_pga_ref) max(MSSM_pga_ref)];

pga_ratio_tmp=reshape(pga_ratio{pp},[length(malawi_lon_site_vec)*length(malawi_lat_site_vec),1]);
pga_ratio_analysis(pp,1:3) =[median(pga_ratio_tmp) mean(pga_ratio_tmp) max(pga_ratio_tmp)];

pga_ratio_tmp_ref=reshape(pga_ratio_ref{pp},[length(malawi_lon_site_vec)*length(malawi_lat_site_vec),1]);
pga_ratio_analysis_ref(pp,1:3) =[median(pga_ratio_tmp_ref) mean(pga_ratio_tmp_ref) max(pga_ratio_tmp_ref)];


 end
 
%% Plot histograms for the spread of hazard map values differences 
 
pp=1; %chose PoE level
 
figure(1000);
 
tiledlayout(2,1,'tilespacing','tight');
nexttile %plot differences between MSSM and P2017 maps 
 
histogram(tmp2{pp},'BinWidth' ,0.025);
xlim([-0.1 0.4]); ylim([0 400]); 
xlabel('Difference in seismic hazard value (g)'); ylabel('Count');

  t4=title('(a)','Fontsize',13,'fontweight','normal');
   set(t4, 'horizontalAlignment', 'left');
   set(t4, 'units', 'normalized');
   h1 = get(t4, 'position');
   set(t4, 'position', [-0.1 h1(2) h1(3)]);

nexttile %plot differences between MSSM and H2015 maps 
 
histogram(h15_tmp2{pp},'BinWidth' ,0.025);
xlim([-0.1 0.4]); ylim([0 400]);
xlabel('Difference in seismic hazard value (g)'); ylabel('Count');

  t4=title('(b)','Fontsize',13,'fontweight','normal');
   set(t4, 'horizontalAlignment', 'left');
   set(t4, 'units', 'normalized');
   h1 = get(t4, 'position');
   set(t4, 'position', [-0.1 h1(2) h1(3)]);



set(gcf,'Position',[440 153 516 644]);

