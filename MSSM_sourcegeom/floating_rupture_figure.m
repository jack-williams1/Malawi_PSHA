% RANDOMLY FLOAT RUPTURES ON FAULT PLANE

% Ruptures from Adapted MSSM event catalog 
% Rupture position constrain by fault length and section width
% Test version of source2site_v3
% Runs for all events in a for loop, and can be used to plot mo_rate
% distribution around fault 
% Used as Fig 5 in PSHA paper

clear all
%load MSSM_Catalog_Adapted running for mutliple faults

mydir  = pwd; idcs   = strfind(mydir,'/');
addpath(mydir(1:idcs(end)-1));
addpath([mydir(1:idcs(end)-1) '/syncat_adaptedMSSM']);
addpath([mydir(1:idcs(end)-1) '/misc_functions']);

load MSSM_sources
load fault_geom_1
load fault_geom_2
load MSSM_catalog_adapted_em


% MSSM source being plotted
source_id = 301; %ID of single source being assessed
event_num1=0; event_num2=0;

%Read in MAFD shape file, update if necessary
MAFD = shaperead('MAFD.shp');
num_traces = length(MAFD);

ellipsoid = almanac('earth','wgs84','meters');
utms = defaultm('utm'); utms.zone = '36L'; utms.geoid = ellipsoid; 
utms.flatlimit = []; utms.maplatlimit = []; utms = defaultm(utms);

%find source(s) in the MAFD
id_xref=readtable('id_xref.xlsx');
if source_id<600
    tmp=find(source_id==id_xref.unofficial_id);
    MAFD_id=id_xref.MAFD_id(tmp);
else %multifault source
    tmp1=find(source_id==num_fault(:,23));
    for gg=1:length(tmp1)
        tmp2=find(num_fault(tmp1(gg),1)==id_xref.unofficial_id);
        MAFD_id(gg)=id_xref.MAFD_id(tmp2);
    end
end

%% Associate ruptures with grid points from faultgrid_wgs84

%Loop through each StochasticEventCatalog
for ii=1:length(StochasticEventCatalog)
    
     gp_tmp={};
     
    %Find events associated with source_id and M>5.4 in event catalogs
    source_indx=find(StochasticEventCatalog{ii}(:,6)==source_id);
    
    for ijk=1:length(source_indx)
        
        syncat_ijk=StochasticEventCatalog{ii}(source_indx(ijk),[1,2,4,6,5]);
        
        if syncat_ijk(3)<5.4 %treat small events as point sources
            [gp]=source2site_v2(syncat_ijk,faultgrid_wgs84_1,faultgrid_wgs84_2,grid_indx1,grid_indx2);
        else
            [gp]=source2site_v3(syncat_ijk,faultgrid_wgs84_1,faultgrid_wgs84_2,num_fault,num_multi_fault,num_MSSM,grid_indx1,grid_indx2);
        end
        gp_tmp{ijk}=gp;
        clear gp
    end
    
    plot_gp{ii}=gp_tmp;
    
    if ii==1 || ii==2
    event_num1=event_num1+length(plot_gp{ii});
    else
    event_num2=event_num2+length(plot_gp{ii});
    end
    
    clear source_indx gp_tmp
end

%% Find number of events each gp of fault surface has hosted

flt_indx = find(num_MSSM(:,1)==source_id); 

if source_id>600
    flt_all_indx = find(source_id==num_fault(:,23),1);  
else
    flt_all_indx = find(num_fault(:,1)==source_id);%for indexing num_fault
end
dip=num_fault(flt_all_indx,8);

grd_indx_f1 = find(grid_indx1==source_id); grd_indx_f2 = find(grid_indx2==source_id); 
source_pts1 = single(faultgrid_wgs84_1{grd_indx_f1});
source_pts2 = single(faultgrid_wgs84_2{grd_indx_f2});

count_1=zeros(length(source_pts1),1); count_2=zeros(length(source_pts2),1);

tmp_count =1;
clock_old = cputime;

for ii=1:length(StochasticEventCatalog)
    
    % source_pts model depends on fault width assumption
    if ii==1 || ii==2
        source_pts=source_pts1;
    else
        source_pts=source_pts2;
    end
    
     for ijk=1:length(plot_gp{ii})
         
        gp=cell2mat(plot_gp{ii}(ijk));

        for ff=1:height(gp)
        
            tmp=find(single(gp(ff,1:3))==source_pts(:,1:3),1);
        
            if (isempty(tmp)==0 && ii==1) || (isempty(tmp)==0 && ii==2)
            count_1(tmp)=count_1(tmp)+1;%Count event for grid point
        
            elseif (isempty(tmp)==0 && ii==3) || (isempty(tmp)==0 && ii==4)
            count_2(tmp)=count_2(tmp)+1;%Count event for grid point
            clear tmp 
            end %end if statement
        
        end %end ff loop
   clear gp
    end %end ijk loop for each event in catalog
    
end %end ii loop for each catalog


%% Create grid of fault in x*y space


f_id = unique(source_pts1(:,4));  mssm_indx =  find(num_MSSM(:,1)==source_id);%for indexing MSSM

%index fault or multifault parts in fault geometry points

for i=1:length(f_id)
    for j=1:length(fault_geometry_points_1)
        tmp=find(f_id(i)==fault_geometry_points_1{j}(1,4));
        if isempty(tmp)==0
           geom_indx(i)=j;
           break
        end
        
    end
end

source_pts1=double(source_pts1); source_pts2=double(source_pts2);
cmap = crameri('batlow');

z1_1={}; z1_2={}; z2_1={}; z2_2={};

for ff=1:length(f_id)
    f_indx1 = find(source_pts1(:,4)==f_id(ff));
    f_indx2 = find(source_pts2(:,4)==f_id(ff));
    [xq_1{ff}, yq_1{ff}]=meshgrid(min(source_pts1(f_indx1,2)):0.002:max(source_pts1(f_indx1,2)),min(source_pts1(f_indx1,1)):0.002:max(source_pts1(f_indx1,1)));
    [xq_2{ff}, yq_2{ff}]=meshgrid(min(source_pts2(f_indx2,2)):0.002:max(source_pts2(f_indx2,2)),min(source_pts2(f_indx2,1)):0.002:max(source_pts2(f_indx2,1)));
    
    z1_1{ff} =  griddata(source_pts1(f_indx1,2), source_pts1(f_indx1,1), count_1(f_indx1,1), xq_1{ff}, yq_1{ff});%Assign grid moment_rate 
    z1_2{ff}= griddata(source_pts1(f_indx1,2), source_pts1(f_indx1,1), source_pts1(f_indx1,3), xq_1{ff}, yq_1{ff}); %Assign grid depth values
    z2_1{ff} =  griddata(source_pts2(f_indx2,2), source_pts2(f_indx2,1), count_2(f_indx2,1), xq_2{ff}, yq_2{ff});
    z2_2{ff}= griddata(source_pts2(f_indx2,2), source_pts2(f_indx2,1), source_pts2(f_indx2,3), xq_2{ff}, yq_2{ff}); 
    clear f_indx1 f_indx2
end

figure(505);

subplot(1,2,1)
for ff=1:length(f_id)
surf(xq_1{ff}, yq_1{ff}, z1_2{ff}, z1_1{ff},'EdgeColor','None'); colormap(cmap); caxis([0 max(count_2(:,1))]); hold on;
end

%cbar=colorbar; cbar.Label.String = 'Number of events'; 
%string_1 = {['Total events: ' num2str(event_num1)]}; hold on
%set(cbar.Title,{'String','Rotation','Position'},{string_1,0,[0 300]}); hold on;
%set(cbar,'position',[.48 .35 .015 .3]); hold on
for g=1:length(geom_indx)
plot3(fault_geometry_points_1{geom_indx(g)}(:,2),fault_geometry_points_1{geom_indx(g)}(:,1), fault_geometry_points_1{geom_indx(g)}(:,3),'k-'); hold on
plot3(fault_geometry_points_2{geom_indx(g)}(:,2),fault_geometry_points_2{geom_indx(g)}(:,1), fault_geometry_points_2{geom_indx(g)}(:,3),'k-'); hold on
end

%Plot MAFD trace of fault
for kk = 1:length(MAFD_id)    
    [LAT,LON] = minvtran(utms,MAFD(MAFD_id(kk)).X,MAFD(MAFD_id(kk)).Y);
    plot3(LON,LAT,zeros(length(LON),1),'r-','linewidth',2); hold on; 
end

axis equal; daspect([0.01 0.01 1.4]); zlim([0 35]); 
xlabel('Longitude'); ylabel('Latitude'); zlabel('Depth (km)');hold on,
ax=gca; ax.ZDir = 'reverse'; ax.XDir = 'reverse';   set(get(gca,'xlabel'),'rotation',+110); 
set(gca,'fontsize',13); view((num_MSSM(mssm_indx,4)-125),(95-dip)) %View parallel to fault strike

subplot(1,2,2)
for ff=1:length(f_id)
surf(xq_2{ff}, yq_2{ff}, z2_2{ff}, z2_1{ff},'EdgeColor','None'); colormap(cmap); hold on% caxis([0 max(mo_r2(:,1))]); hold on;
end

cbar=colorbar; cbar.Label.String = 'Number of events'; set(cbar,'position',[.935 .35 .015 .3]); hold on
string_2 = {['Length-limited {\it W} total events: ' num2str(event_num1) newline 'Layer-limited {\it W} total events: ' num2str(event_num2)]}; hold on
set(cbar.Title,{'String','Rotation','Position'},{string_2,0,[-13 300]}); hold on;
for g=1:length(geom_indx)
plot3(fault_geometry_points_2{geom_indx(g)}(:,2),fault_geometry_points_2{geom_indx(g)}(:,1), fault_geometry_points_2{geom_indx(g)}(:,3),'k-'); hold on
end

for kk = 1:length(MAFD_id)    
    [LAT,LON] = minvtran(utms,MAFD(MAFD_id(kk)).X,MAFD(MAFD_id(kk)).Y);
    plot3(LON,LAT,zeros(length(LON),1),'r-','linewidth',2); hold on; 
end

axis equal; daspect([0.01 0.01 1.4]); xlabel('Longitude'); ylabel('Latitude'); zlabel('Depth (km)');hold on,
ax=gca; ax.ZDir = 'reverse'; ax.XDir = 'reverse';  set(get(gca,'xlabel'),'rotation',+110);
set(gca,'fontsize',13); view((num_MSSM(mssm_indx,4)-125),(95-dip))

set(gcf, 'Position', [658 122 1351 863])

