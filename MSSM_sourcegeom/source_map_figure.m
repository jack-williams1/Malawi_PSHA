%% Plot figure of MSSM sources stored in fault_geom

% Plot fault geom_1 (1) or fault geom2 (2)
geom_option =2;

if geom_option ==1
    load fault_geom_1 
    faultgrid_wgs84=faultgrid_wgs84_1;
    fault_geometry_points=fault_geometry_points_1;
    grid_indx=grid_indx1;
else
    load fault_geom_2
    faultgrid_wgs84=faultgrid_wgs84_2;
    fault_geometry_points=fault_geometry_points_2;
    grid_indx=grid_indx2;
end

%% Load functions needed for fault geometry calculations
mydir  = pwd; idcs   = strfind(mydir,'/');
addpath(mydir(1:idcs(end)-1));
addpath([mydir(1:idcs(end)-1) '/misc_functions']);
addpath([mydir(1:idcs(end)-1) '/gis_files']);

load map_data_EastAfrica
load MSSM_sources

%load in Malawi Active Fault Database (MAFD)
MAFD = shaperead('MAFD.shp');

%Name of faults in the MAFD but not in the MSSM (see Table S2 in MSSM manuscript)
removed_faults=["Nchalo","Mudi","Jimbe","Chileka","Nguluwe","Lirangwe River","Linjidzi","Ngondo-1","Ngondo-2","Katundu"...
    "Namiyala-1","Namiyala-2","Namiyala-3","Chilongwelo","Leopard Bay-2","South Basin Fault 4","Central Basin Fault 4",...
    "Central Basin Fault 9","Central Basin Fault 10","Central Basin Fault 22","Hara Plain","South Karonga East","Lupaso"];

num_traces = length(MAFD);

load('syncat_PSHA_MSSM_input','Region')

ellipsoid = almanac('earth','wgs84','meters');
utms = defaultm('utm'); utms.zone = '36L'; utms.geoid = ellipsoid; 
utms.flatlimit = []; utms.maplatlimit = []; utms = defaultm(utms);

s_th = 35; % Thickness of seismogenic crust


%% Remove duplicate faults 

%Faults to remove

%Thyolo 312
%Wamkurumadzi 321
%Lisungwe 319
%Chirobwe-Ncheu 327
%Bilia-Mtakataka 329
%Metangula 337
%SB Fault 2b 348
%SB Fault 5b 353
%SB Fault 5c 354
%SB Fault 7a 359
%SB Fault 11b 350
%SB Fault 12a 361
%Sabi-2 396
%South Karonga-2 395
%Mbiri-2 399
%North Basin 4a 416
%North Basin 9b 413
%North Basin 15a 417 
dupfault=[312 319 321 327 329 337 348 353 354 359 350 361 396 395 399 413 416 417];

for i=1:length(faultgrid_wgs84)
    if ismember(faultgrid_wgs84{i}(1,4),dupfault)==1
        faultgrid_wgs84{i}=[0 0 0 0 0 0 0];
    end
end

tmp=faultgrid_wgs84';
faultgrid_wgs84_num=cell2mat(tmp);

%% Randomly assign color to each grid point based on section number
%Use jet color scheme up to value of fault with most sections
color=parula(max(faultgrid_wgs84_num(:,6)));
tmp_a=[1:1:length(color)];

%Remove multi faults that have been combined
tmp_idx=find(faultgrid_wgs84_num(:,1)==0);
faultgrid_wgs84_num(tmp_idx,:)=[];

% Fix random number, only needs to run once
% Ensures fault sections always given same colors
%fix_rng=rng; %save('fix_rng','fix_rng')

load fix_rng; rng=fix_rng;

for i=1:length(faultgrid_wgs84)
     %Ignore faults who have been combined into multi_faults
    if faultgrid_wgs84{i}(1,2)>0
        
    %Generate random intergers for each fault, whose length is equal 
    %to the number of sections
    col{i}(:,2)=(randperm(length(unique(faultgrid_wgs84_num(:,6))),length(unique(faultgrid_wgs84{i}(:,5)))))';
  %  col{i}(:,1)=[1:1:length(unique(faultgrid_wgs84{i}(:,5)))];
    col{i}(:,1)=unique(faultgrid_wgs84{i}(:,5));
        for gg=1:length(unique(faultgrid_wgs84{i}(:,5)))
            for ff=1:length(faultgrid_wgs84{i})
    %For each fault grid point, find the random number based on its index            
               if faultgrid_wgs84{i}(ff,5) == col{i}(gg,1)
                    grid_color_index{i}(ff)=col{i}(gg,2);
                end
        end  
    end   
    end
end

grid_color_index_num=cell2mat(grid_color_index);
grid_color_index_num = grid_color_index_num';

%Assign the grid points its color based on its random index
grid_color=zeros(length(faultgrid_wgs84_num),5);
for i=1:length(faultgrid_wgs84_num)
        
    for ll=1:length(color)
        if grid_color_index_num(i)==tmp_a(ll);
             grid_color(i,2:4)=color(ll,:); 
             grid_color(i,1)=faultgrid_wgs84_num(i,6);
        end
    end
end


%% Plot fault grids

LakeMalawi = shaperead('malawi_lake.shp');
LakeMalawiCoord = [LakeMalawi.Y(1,1:end-1)' LakeMalawi.X(1,1:end-1)'];

figure(234);

%Plot grids with randomly assigned color
% Need to find a way to make this quicker and not run on a for loop

simu_count=1; clock_old = cputime;

for i=1:length(grid_color)
    
     if i== 10000*simu_count
        clock_new = cputime;
        disp(['Simulation cycle is: ',num2str(i),' out of ',num2str(length(grid_color)),' & Required time is: ',num2str(clock_new-clock_old)]);
        clock_old = clock_new;
        simu_count = simu_count + 1;
     end
    
plot3(faultgrid_wgs84_num(i,2),faultgrid_wgs84_num(i,1),-1*faultgrid_wgs84_num(i,3),'.','color',grid_color(i,2:4)) 
hold on
end

hold on


%Plot MAFD traces
for kk = 1:num_traces 
    [LAT,LON] = minvtran(utms,MAFD(kk).X,MAFD(kk).Y);
    if contains(MAFD(kk).fault_name,removed_faults)==1
        %map faults not included in the MSSM as black
        plot3(LON,LAT,ones(length(LON),1),'k-','linewidth',2); hold on;
    else
        plot3(LON,LAT,ones(length(LON),1),'r-','linewidth',1.5); hold on; 
    end
end

ylabel('Latitude','FontSize',16); xlabel('Longitude','FontSize',16); zlabel('Depth (km)');
set(get(gca,'ylabel'),'rotation',50); set(get(gca,'xlabel'),'rotation',335);zticks([-35 0])
zticklabels({'-35','0'});
%Plot DEM, border & lake

%hillshade generated from GDEM in gis_file folder with stored xy data
malawi_hillshade = imread('malawi_hillshade.jpg'); 
load hillshade_xy_data %load coordinate locations
floor_image = hgtransform('Matrix',makehgtform('translate',[0 0 -35]));
image('XData',Lon_Malawi,'YData',Lat_Malawi,'CData',malawi_hillshade,'Parent',floor_image); colormap(gray); 
plot3(MapData2(:,2),MapData2(:,1),-34.9*ones(length(MapData2(:,1)),1),'w-') ; hold on  
fill3(LakeMalawiCoord(:,2),LakeMalawiCoord(:,1),-34.9*ones(length(LakeMalawiCoord(:,1)),1),'w') ; hold on  

xlim([Region(1) Region(2)]); ylim([Region(3) Region(4)]); grid on;
ax=gca; daspect([0.01 0.01 1]); zlim([-35 1]); hold on

%view(0,90) %birdseye view
view(30,45) %oblique view
set(gca,'fontsize',13)
set(gcf,'position',[681 151 980 828])
