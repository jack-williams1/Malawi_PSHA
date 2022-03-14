%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%   Probabilistic Seismic Hazard Analysis for Malawi   %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tic
% Seismic hazard map for a single ground motion parameter
% FOR ALL PSHA CATALOGS, RUN SEPERATELY ONLY
% Splits catalogs to run on BlueCrystal

clear
close all

clock_start = cputime;

num_par=3; %script number in fake parallelisation

mydir  = pwd; idcs   = strfind(mydir,'/');
addpath(mydir(1:idcs(end)-1));

load('syncat_PSHA_MSSD_input','num_GMPEindex','GMPEindex','weight_GMPEindex','FMtype','FMweight','Z10_Z25','vs30_site_ref',...
    'obs_duration','lower_return_period','num_need','max_bg','par_opts','syncat_name',...
    'Region','PSHA_Zone','site_corner_SW','site_corner_NE','site_grid_interval','T_map');

addpath([mydir(1:idcs(end)-1) '/syncat_MSSD']); addpath([mydir(1:idcs(end)-1) '/syncat_adaptedMSSD']); 
addpath([mydir(1:idcs(end)-1) '/syncat_bg']); addpath([mydir(1:idcs(end)-1) '/misc_functions']);
addpath([mydir(1:idcs(end)-1) '/gis_files']); addpath([mydir(1:idcs(end)-1) '/MSSD_sourcegeom']);
addpath([mydir(1:idcs(end)-1) '/GMPE']);
load MSSD_sources

vs30_data = textread('malawi_Vs30_active.txt'); 

lat_site_vec = site_corner_SW(1):site_grid_interval:site_corner_NE(1);
lon_site_vec = site_corner_SW(2):site_grid_interval:site_corner_NE(2);

[Lat_site,Lon_site] = meshgrid(lat_site_vec,lon_site_vec);

site = [Lat_site(:) Lon_site(:)];

num_site = size(site,1);

%% Sort vs30 data

vs30_site =zeros(length(num_site),2);

for ii = 1:num_site
    [~,mindiff_id] = min((vs30_data(:,2)-site(ii,1)).^2 + (vs30_data(:,1)-site(ii,2)).^2);
    %Lat and Long increments for vs30 data different for site, so find the vs30
    %data that is closest to each site increment
    vs30_site(ii,1) = vs30_data(mindiff_id,3);
end

vs30_site(:,2) = ones(num_site,1)*vs30_site_ref;


%% GMPE options 
% 1)  Boore & Atkinson (2008) for shallow crustal events (NGA-WEST1)
% 2)  Campbell & Bozorgnia (2008) for shallow crustal events (NGA-WEST1)
% 3)  Chiou & Youngs (2008) for shallow crustal events (NGA-WEST1)
% 4)  Akkar & Bommer (2010) - Mw = 5.0-7.6; constant sigma model
% 5)  Boore, Stewart, Seyhan & Atkinson (2014) for shallow crustal events (NGA-WEST2)
% 6)  Akkar, Sandikkaya, & Bommer (2014) for shallow crustal events (Rjb model)
% 7)  Akkar, Sandikkaya, & Bommer (2014) for shallow crustal events (Repi model)
% 8)  Akkar, Sandikkaya, & Bommer (2014) for shallow crustal events (Rhypo model)
% 9)  Cauzzi, Faccioli, Vanini, & Bianchini (2015) for shallow crustal events
% 10) Chiou & Youngs (2014) for shallow crustal events
% 11) Atkinson and Adams (2013) for stable continental events (6th Seismic Hazard Model of Canada) - best/lower/upper branches 
% 12) Goulet et al. (2017) for stable continental events (NGA-East preliminary with site amplification; 6th Seismic Hazard Model of Canada) - models 1 to 13  

% Spectral vibration periods; "T = 0" indicates the PGA.

for ii = 1:num_GMPEindex
    COEF = GMPEcoef_Malawi(GMPEindex(ii),T_map);
    assignin('base',strcat('COEF',num2str(ii)),COEF);
    clear COEF
end

%% Probabilistic seismic hazard analysis for background sources
% All background earthquakes are treated as point sources

disp('                                       ');
disp('- - - PSHA for background sources - - -');
 
load syncat_bg

%BG_catalog
%1) Dist from Malawi Region (+ve =outside region)
%2) Latitude 
%3) Longitude
%4) Magnitude
%5) Time

%Remove M>7 bg events from within Malawi
rem_event=find(bg_catalog(:,4) >max_bg & bg_catalog(:,1)<0);
bg_catalog(rem_event,:)=[];


%COMMENT OUT BEFORE RUNNING CODES PROPERLY ON BLUE CRYSTAL
%bg_catalog=bg_catalog(1:round(length(bg_catalog)/10),:);

global_count_bg = length(bg_catalog);

%Divide bg_catalog for parallelisation
divide_bg_indx = [floor(global_count_bg*((par_opts(num_par)-1)/length(par_opts)))+1 floor(global_count_bg*(par_opts(num_par)/length(par_opts)))];
bg_catalog_divide=bg_catalog([divide_bg_indx(1):divide_bg_indx(2)],:);

tmp_global_count_bg=length(bg_catalog_divide);
bg_catalog_indx=(divide_bg_indx(1):1:divide_bg_indx(2))';

%Focal depth parameters
%1 Mean value
%2 Standard deviation
%3 Min value
%4 Max value

fd_para = [20 5 5 35;]; 
fd = min(max(fd_para(1)+fd_para(2)*randn(tmp_global_count_bg,1),fd_para(3)),fd_para(4));

ran = rand(tmp_global_count_bg,1);

dip_bg = [53 90]; bg_fm=zeros(length(tmp_global_count_bg),2);
ran_fm =  rand(tmp_global_count_bg,1);

%assign events to be randomly normal or strike-slip
for i=1:2
    tmp=find(ran > sum(FMweight(1:i-1))  & ran <= sum(FMweight(1:i)));
    bg_fm(tmp,1)=FMtype(i);
    bg_fm(tmp,2)=dip_bg(i);
end

clock_old = cputime;
GM_bg_em = cell(num_GMPEindex,1); GM_bg_em_ref = cell(num_GMPEindex,1);

%% Perform PSHA

for jj = 1:num_GMPEindex
    
    group_count = 1; count=0; count_gm_bg = 0; tmp_count_bg = 1; 
    
    GM_bg = single(zeros(num_need,num_site)); GM_bg_ref = single(zeros(num_need,num_site));
    max_GM_bg =single(zeros(num_need,num_site)); max_ID_bg =single(zeros(num_need,num_site));
    max_GM_bg_ref =single(zeros(num_need,num_site)); max_ID_bg_ref =single(zeros(num_need,num_site));    
    
    assignin('base','COEF',eval(strcat('COEF',num2str(jj))));
   
    
    for ii = 1:tmp_global_count_bg %loops through all events in catalog
        
        count_gm_bg = count_gm_bg + 1; count=count+1;
        
        if count_gm_bg == tmp_count_bg*10000
            tmp_count_bg = tmp_count_bg + 1;
            disp(['Current bg event number is : ',num2str(count_gm_bg),' out of ',num2str(tmp_global_count_bg), ' for GMPE ' num2str(jj)]);
        end
        
        repi  = deg2km(distance(site(:,1),site(:,2),bg_catalog_divide(ii,3), bg_catalog_divide(ii,2))); %Determines source to site distance for each event for each site
        rhypo = sqrt(repi.^2 + fd(ii)^2); 
        
        %COEF is structure for GMPEindex 11 & 12 and so must be assessed within function
          if GMPEindex(jj)== 11 || GMPEindex(jj)== 12
        [GMmedtmp,Sintratmp,Sintertmp] = GMPE_Malawi(GMPEindex(jj),COEF,bg_catalog_divide(ii,4),repi,rhypo,0,repi,rhypo,vs30_site(:,1),bg_fm(ii,1),fd(ii),Z10_Z25,bg_fm(ii,2),T_map); 
        %For ref vs30 value
        [GMmedtmp_ref,Sintratmp_ref,Sintertmp_ref] = GMPE_Malawi(GMPEindex(jj),COEF,bg_catalog_divide(ii,4),repi,rhypo,0,repi,rhypo,vs30_site(:,2),bg_fm(ii,1),fd(ii),Z10_Z25,bg_fm(ii,2),T_map);   
        
          else
        %For USGS vs30 value
        [GMmedtmp,Sintratmp,Sintertmp] = GMPE_Malawi(GMPEindex(jj),ones(num_site,1)*COEF,bg_catalog_divide(ii,4),repi,rhypo,0,repi,rhypo,vs30_site(:,1),bg_fm(ii,1),fd(ii),Z10_Z25,bg_fm(ii,2),T_map); 
        %For ref vs30 value
        [GMmedtmp_ref,Sintratmp_ref,Sintertmp_ref] = GMPE_Malawi(GMPEindex(jj),ones(num_site,1)*COEF,bg_catalog_divide(ii,4),repi,rhypo,0,repi,rhypo,vs30_site(:,2),bg_fm(ii,1),fd(ii),Z10_Z25,bg_fm(ii,2),T_map);
          end
          
        GM_bg(count,:) = GMmedtmp.*exp(sqrt(Sintratmp.^2+Sintertmp.^2).*(randn(num_site,1)));%The ground motion at each site, from each event
        GM_bg_ref(count,:) = GMmedtmp_ref.*exp(sqrt(Sintratmp_ref.^2+Sintertmp_ref.^2).*(randn(num_site,1)));%The ground motion at each site, from each event
        
        if ii == group_count*num_need || ii==tmp_global_count_bg%Number of events in group size
            
             for kk = 1:num_site
    
            tmp_GMstore = sortrows(vertcat([GM_bg(1:count,kk) bg_catalog_indx((group_count-1)*num_need+1:min(group_count*num_need,tmp_global_count_bg))], [max_GM_bg(:,kk) max_ID_bg(:,kk)]),1);
            tmp_GMstore_ref = sortrows(vertcat([GM_bg_ref(1:count,kk) bg_catalog_indx((group_count-1)*num_need+1:min(group_count*num_need,tmp_global_count_bg))],[max_GM_bg_ref(:,kk) max_ID_bg_ref(:,kk)]),1);
    
            %store num_need largest ground motions for each site
            max_GM_bg(1:num_need,kk) = tmp_GMstore(length(tmp_GMstore )-num_need+1:length(tmp_GMstore),1);
            %max_ID_bg(1:num_need,kk) = tmp_GMstore(length(tmp_GMstore)-num_need+1:length(tmp_GMstore),2);
            max_GM_bg_ref(1:num_need,kk) = tmp_GMstore_ref(length(tmp_GMstore_ref)-num_need+1:length(tmp_GMstore_ref),1);
            %max_ID_bg_ref(1:num_need,kk) = tmp_GMstore_ref(length(tmp_GMstore_ref)-num_need+1:length(tmp_GMstore_ref),2);
    
             end %end for loop for num_site
             
        group_count = group_count + 1; count=0;
        GM_bg = single(zeros(num_need,num_site)); GM_bg_ref = single(zeros(num_need,num_site));      
        end %end if statement for group_size==ii
          
    end %end ii loop
  
    %Store GM for each GMPE
    %Note event IDs not currently stored, but this can be added back in
    GM_bg_em{jj} = max_GM_bg; GM_bg_em_ref{jj} = max_GM_bg_ref;
    
end

%% Save parameter based on num_para
%Note stored as single variables to save memory
assignin('base',strcat('AMAX_GM_bg_em_',string(num_par)),GM_bg_em)
assignin('base',strcat('AMAX_GM_bg_em_ref_',string(num_par)),GM_bg_em_ref)

save(strcat('AMAX_GM_bg_',string(num_par)),strcat('AMAX_GM_bg_em_',string(num_par)),strcat('AMAX_GM_bg_em_ref_',string(num_par)),'-v7.3');

%% Load fault based catalog and geometry

%load fault geometrical models for sourc2site calc
load fault_geom_1; load fault_geom_2; load sec_geom_1;

%load fault catalogs
load(char(syncat_name),'EQCAT');
load('MSSD_Catalog_Adapted_em','StochasticEventCatalog');

%from all MSSD direct (EQCAT) and/or MSSD adapted take
 
 % 1) event number
 % 2) occurrence time
 % 3) magnitude
 % 4) source_id
 % 5) width case (MSSD adapted only)
%and combine into one cell array

syncat{1}=EQCAT(:,[1,2,5,7]);%MSSD direct catalog
%COMMENT FOR TESTING, REMOVE WHEN ON BLUE CRYSTAL
%syncat{1}=syncat{1}(1:round(length(syncat{1})/10),:);

global_count=length(syncat{1});

divide_indx = [floor(global_count*((par_opts(num_par)-1)/length(par_opts)))+1 floor(global_count*(par_opts(num_par)/length(par_opts)))];
syncat{1}=syncat{1}([divide_indx(1):divide_indx(2)],:);
syncat_indx{1}=(divide_indx(1):1:divide_indx(2))';

for i=1:length(StochasticEventCatalog)
    syncat{i+1}=StochasticEventCatalog{i}(:,[1,2,4,6,5]);%update with more event catalogs as they are generatred
    
    %COMMENT FOR TESTING, REMOVE WHEN ON BLUE CRYSTAL
    %syncat{i+1}=syncat{i+1}(1:round(length(syncat{i+1})/10),:);
    
    global_count=length(syncat{i+1});
    divide_indx = [floor(global_count*((par_opts(num_par)-1)/length(par_opts)))+1 floor(global_count*(par_opts(num_par)/length(par_opts)))];
    syncat{i+1}=syncat{i+1}([divide_indx(1):divide_indx(2)],:);
    syncat_indx{i+1}=(divide_indx(1):1:divide_indx(2))';
    
end


%% Run fault-based PSHA

for ss=1:length(syncat)

    
global_count=length(syncat{ss});
max_GM_fs =single(zeros(num_need,num_site,num_GMPEindex)); max_ID_fs =single(zeros(num_need,num_site,num_GMPEindex));
max_GM_fs_ref =single(zeros(num_need,num_site,num_GMPEindex)); max_ID_fs_ref =single(zeros(num_need,num_site,num_GMPEindex));
ran = rand(global_count,1); 
clock_old = cputime;

group_count = 1; count=0; tmp_count = 1; 
    
GM_fs = single(zeros(num_need,num_site,num_GMPEindex)); GM_fs_ref = single(zeros(num_need,num_site,num_GMPEindex));

for ijk = 1:global_count

   if ijk == tmp_count*5000
      tmp_count = tmp_count + 1;
      clock_new = cputime;
      disp(['Current event number is : ',num2str(ijk),' out of ',num2str(global_count),...
      ' and syncat ', num2str(ss), ' & Required time is: ',num2str(clock_new-clock_old)]);
      clock_old = clock_new;
   end
   
   count=count+1;
   syncat_ijk=syncat{ss}(ijk,:);
   
    %directly obtain grids from known earthquake geometry for direct MSSD
    if ss==1
        [gp]=source2site_v1(syncat_ijk,faultgrid_wgs84_1,secgrid_wgs84_1,num_sec,num_fault,grid_indx1);
    else
        
        %float gp on fault grids depending on magnitude for adapted MSSD
        %treat small events as point sources
        if syncat_ijk(3)<5.4
            [gp]=source2site_v2(syncat_ijk,faultgrid_wgs84_1,faultgrid_wgs84_2,grid_indx1,grid_indx2);
            
        %treat large events as 2D planes  
        else
            [gp]=source2site_v3(syncat_ijk,faultgrid_wgs84_1,faultgrid_wgs84_2,num_fault,num_multi_fault,num_MSSD,grid_indx1,grid_indx2);
        end
    
    end
    
    %index source dip from MSSD_sources
    if syncat_ijk(4)>600
            %find dip angle from num_multi_fault from first faults
            dip = num_fault(find(syncat_ijk(4)==num_fault(:,23),1,'first'),8);      
    elseif syncat_ijk(4)>300
             %find dip angle from num_fault for single section faults
             dip = num_fault(find(syncat_ijk(4)==num_fault(:,1)),8);       
    else
             %find dip angle from num_fault for single section ruptures    
             dip = num_sec(find(syncat_ijk(4)==num_sec(:,1)),10);
    end
    
    if size(gp,1)==1 %if event treated as point source
        
    repi = deg2km(distance(site(:,1),site(:,2),gp(1),gp(2))); % Rjb, 2D distance
    rhypo = sqrt(repi.^2 + gp(:,3).^2); % Rcd, 3D distance   
   
    DIST = [repi rhypo repi rhypo];
    
    else %for floating ruptures, have to iteratively run through each site

    DIST=zeros(num_site,4);
    
        for ii=1:num_site
            %Returns four types of distances
            %DIST(:,1) = Rjb, 2D distance
            %DIST(:,2) = Rcd, 3D distance
            %DIST(:,3) = Repi, epicentral distance (random point on top of fault plane)
            %DIST(:,4) = Rhypo, hypocentral distance (random point in fault plane) 
            dist_fault1 = deg2km(distance(site(ii,1),site(ii,2),gp(:,1),gp(:,2))); % Rjb, 2D distance
            dist_fault2 = sqrt(dist_fault1.^2 + gp(:,3).^2); % Rcd, 3D distance
         
            ran_dist   = ceil(size(gp,1)*rand(1,1));
            repi  = dist_fault1(ran_dist); %epicentral distance (random point on top of fault plane)
            rhypo = dist_fault2(ran_dist); %hypocentral distance (random point in fault plane) 
            DIST(ii,:)= [min(dist_fault1); min(dist_fault2); repi; rhypo]; 
        end
        
    end %end if statement for point source
    
    %for simulated event, loop through each GMPE
    for jj = 1:num_GMPEindex
        
    assignin('base','COEF',eval(strcat('COEF',num2str(jj))));
   
    %COEF is structure for GMPEindex 11 & 12 and so must be assessed within function
      if GMPEindex(jj)== 11 || GMPEindex(jj)== 12
    %For USGS vs30 value
        [GMmedtmp,Sintratmp,Sintertmp] = GMPE_Malawi(GMPEindex(jj),COEF,syncat_ijk(3),DIST(:,1),DIST(:,2),0,DIST(:,3),DIST(:,4),vs30_site(:,1),FMtype(1),min(gp(:,3)),Z10_Z25,dip,T_map); 
    %For ref vs30 value
        [GMmedtmp_ref,Sintratmp_ref,Sintertmp_ref] = GMPE_Malawi(GMPEindex(jj),COEF,syncat_ijk(3),DIST(:,1),DIST(:,2),0,DIST(:,3),DIST(:,4),vs30_site(:,2),FMtype(1),min(gp(:,3)),Z10_Z25,dip,T_map);   
      else
    %For USGS vs30 value
        [GMmedtmp,Sintratmp,Sintertmp] = GMPE_Malawi(GMPEindex(jj),ones(num_site,1)*COEF,syncat_ijk(3),DIST(:,1),DIST(:,2),0,DIST(:,3),DIST(:,4),vs30_site(:,1),FMtype(1),min(gp(:,3)),Z10_Z25,dip,T_map); 
    %For ref vs30 value
        [GMmedtmp_ref,Sintratmp_ref,Sintertmp_ref] = GMPE_Malawi(GMPEindex(jj),ones(num_site,1)*COEF,syncat_ijk(3),DIST(:,1),DIST(:,2),0,DIST(:,3),DIST(:,4),vs30_site(:,2),FMtype(1),min(gp(:,3)),Z10_Z25,dip,T_map); 
      end
      
        GM_fs(count,:,jj) = GMmedtmp.*exp(sqrt(Sintratmp.^2+Sintertmp.^2).*(randn(num_site,1)));%The ground motion at each site, from each event
        GM_fs_ref(count,:,jj) = GMmedtmp_ref.*exp(sqrt(Sintratmp_ref.^2+Sintertmp_ref.^2).*(randn(num_site,1)));%The ground motion at each site, from each event
    
     if ijk == group_count*num_need || ijk==global_count %Number of events in group size
            
             for kk = 1:num_site
    
            tmp_GMstore = sortrows(vertcat([GM_fs(1:count,kk,jj) syncat_indx{ss}((group_count-1)*num_need+1:min(group_count*num_need,global_count))], [max_GM_fs(:,kk,jj) max_ID_fs(:,kk,jj)]),1);
            tmp_GMstore_ref = sortrows(vertcat([GM_fs_ref(1:count,kk,jj) syncat_indx{ss}((group_count-1)*num_need+1:min(group_count*num_need,global_count))],[max_GM_fs_ref(:,kk,jj) max_ID_fs_ref(:,kk,jj)]),1);
    
            %store num_need largest ground motions for each site
            max_GM_fs(1:num_need,kk,jj) = tmp_GMstore(length(tmp_GMstore)-num_need+1:length(tmp_GMstore),1);
            %max_ID_fs(1:num_need,kk) = tmp_GMstore(length(tmp_GMstore)-num_need+1:length(tmp_GMstore),2);
            max_GM_fs_ref(1:num_need,kk,jj) = tmp_GMstore_ref(length(tmp_GMstore_ref)-num_need+1:length(tmp_GMstore_ref),1);
            %max_ID_fs_ref(1:num_need,kk) = tmp_GMstore_ref(length(tmp_GMstore_ref)-num_need+1:length(tmp_GMstore_ref),2);
    
             end %end for loop for num_site
            
         if jj==num_GMPEindex  %reset variables at end of jj loop  
            group_count = group_count + 1; count=0;
            GM_fs = single(zeros(num_need,num_site,num_GMPEindex)); GM_fs_ref = single(zeros(num_need,num_site,num_GMPEindex)); 
         end
         
      end %end if statement for group_size==ii
            
    end %end jj loop for GMPE
   
end %end ijk loop for syncat  

AMAX_GM_fs{ss}=max_GM_fs; AMAX_ID_fs{ss}=max_ID_fs; 
AMAX_GM_fs_ref{ss}=max_GM_fs_ref; AMAX_ID_fs_ref{ss}=max_ID_fs_ref;

end %end ss loop

%% Save parameter based on num_para
%Note stored as single variables to save memory
assignin('base',strcat('AMAX_GM_fs_',string(num_par)),AMAX_GM_fs)
assignin('base',strcat('AMAX_GM_fs_ref_',string(num_par)),AMAX_GM_fs_ref)

run_time=toc;

save(strcat('AMAX_GM_fs_',string(num_par)),strcat('AMAX_GM_fs_',string(num_par)),strcat('AMAX_GM_fs_ref_',string(num_par)),'-v7.3','run_time')
clock_end = cputime;
disp(['Calculation time: ',num2str(clock_end - clock_start)]);

