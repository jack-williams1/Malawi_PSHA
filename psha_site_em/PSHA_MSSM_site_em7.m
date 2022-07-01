%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%   Probabilistic Seismic Hazard Analysis for Malawi-1   %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tic

% Synthetic earthquake catalog is generated from syncat_MSSM.m
% Epistemic uncertainity assessed through ensemble modelling
% Site specific PSHA, returns ground motions only
% See PSHA_MSSM_site_figures to assess ground motions

clear
close all

%% PSHA set-up

clock_start = cputime;

num_par=7; %script number in fake parallelisation

mydir  = pwd; idcs   = strfind(mydir,'/');
addpath(mydir(1:idcs(end)-1));

addpath([mydir(1:idcs(end)-1) '/syncat_MSSM']); addpath([mydir(1:idcs(end)-1) '/syncat_adaptedMSSM']); 
addpath([mydir(1:idcs(end)-1) '/syncat_bg']); addpath([mydir(1:idcs(end)-1) '/misc_functions']);
addpath([mydir(1:idcs(end)-1) '/gis_files']); addpath([mydir(1:idcs(end)-1) '/MSSM_sourcegeom']);
addpath([mydir(1:idcs(end)-1) '/GMPE']);

%load in PSHA options. To change update excel spreadsheet and then updated
%related Matlab variable
load('syncat_PSHA_MSSM_input','num_GMPEindex','GMPEindex','weight_GMPEindex','Z10_Z25','vs30_site_ref',...
    'obs_duration','lower_return_period','num_need','max_bg','par_opts','FMtype','FMweight',...
    't_limit','NumSimu','T','num_T','site','syncat_name');

load MSSM_sources
vs30_data = textread('malawi_Vs30_active.txt'); 


%% GMPE options:
% 1)  Boore & Atkinson (2008) for shallow crustal events (NGA-WEST1)
% 2)  Campbell & Bozorgnia (2008) for shallow crustal events (NGA-WEST1)
% 3)  Chiou & Youngs (2008) for shallow crustal events (NGA-WEST1)
% 4)  Akkar & Bommer (2010) - Mw = 5.0-7.6; constant sigma model
% 5)  Boore, Stewart, Seyhan & Atkinson (2014) for shallow crustal events (NGA-WEST2)
% 6)  Akkar, Sandikkaya, & Bommer (2014) for shallow crustal events (Rjb model)
% 7)  Akkar, Sandikkaya, & Bommer (2014) for shallow crustal events (Repi model)
% 8)  Akkar, Sandikkaya, & Bommer (2014) for shallow crustal events (Rhypo model)
% 9)  Cauzzi, Faccioli, Vanini, & Bianchini (2015) for shallow crustal events
% 10) Chiou & Youngs (2014) for shallow crustal events (NGA-WEST1)
% 11) Atkinson and Adams (2013) for stable continental events (6th Seismic Hazard Model of Canada) - best/lower/upper branches 
% 12) Goulet et al. (2017) for stable continental events (NGA-East preliminary with site amplification; 6th Seismic Hazard Model of Canada) - models 1 to 13 

%For each GMPE, finds the coefficients for the equation for each spectral
%vibration period
for ii = 1:num_GMPEindex
    COEF = GMPEcoef_Malawi(GMPEindex(ii),T);
    assignin('base',strcat('COEF',num2str(ii)),COEF);% Adds COEF to different objects, COEF1 etc
    clear COEF
end

%% Find site VS30 as input from syncat_PSHA_input

num_site=size(site,1);

%set vs30_site base case to 300 m/s (vs30_site_opp=1) or USGS(vs30_site_op=2)
vs30_site_opp=1;

if vs30_site_opp==1

vs30_site =ones(num_site,1)*300;

else

    vs30_site =zeros(length(num_site),2);
    
    for ii = 1:num_site
        [~,mindiff_id] = min((vs30_data(:,2)-site(ii,1)).^2 + (vs30_data(:,1)-site(ii,2)).^2);
        %Lat and Long increments for vs30 data different for site, so find the vs30
        %data that is closest to each site increment
        vs30_site(ii,1) = vs30_data(mindiff_id,3);
    end
end
    
vs30_site(:,2) = ones(num_site,1)*vs30_site_ref;


%% Probabilistic seismic hazard analysis for background sources
% All background earthquakes are treated as point sources
 
load syncat_bg

%BG_catalog
%1) Dist from Malawi Region +ve =outside region
%2) Latitude
%3) Longitude
%4) Magnitude
%5) Time

%Remove M>7 bg events from within Malawi
rem_event=find(bg_catalog(:,4) >max_bg & bg_catalog(:,1)<0);
bg_catalog(rem_event,:)=[];

%COMMENT/REMOVE OUT BEFORE RUNNING ON BLUE CRYSTAL%
%bg_catalog=bg_catalog(1:500000,:);

global_count_bg=length(bg_catalog);

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

fd = min(max(fd_para(1)+fd_para(2)*randn(tmp_global_count_bg,1),fd_para(3)),fd_para(4));%randomise focal depths for each event between 5-35 km

dip_bg = [53 90]; bg_fm=zeros(length(tmp_global_count_bg),2);
ran_fm =  rand(tmp_global_count_bg,1);

rhypo_bg = zeros(num_site,tmp_global_count_bg); repi_bg = zeros(num_site,tmp_global_count_bg);

%randomly assign events to be normal or strike-slip
for i=1:2
    tmp=find(ran_fm > sum(FMweight(1:i-1))  & ran_fm <= sum(FMweight(1:i)));
    bg_fm(tmp,1)=FMtype(i);
    bg_fm(tmp,2)=dip_bg(i);
end


%% Peform PSHA for each GMPE individually

EPS_bg = single(zeros(tmp_global_count_bg,num_T,num_site));
GM_bg  = single(zeros(num_GMPEindex,tmp_global_count_bg,num_T,num_site));
GM_bg_ref  = single(zeros(num_GMPEindex,tmp_global_count_bg,num_T,num_site)); 
mean_GM_bg = single(zeros(tmp_global_count_bg,num_T,num_site));
mean_GM_bg_ref = single(zeros(tmp_global_count_bg,num_T,num_site));

for gg = 1:num_GMPEindex
        
    assignin('base','COEF',eval(strcat('COEF',num2str(gg))));% Add in GMPE coefficients for selected equation
    
    count_gm_bg = 0;
    tmp_count_bg = 1;

    for ii = 1:tmp_global_count_bg %For each event assigned to a specific GMPE
        
        count_gm_bg = count_gm_bg + 1;
        
        if count_gm_bg == tmp_count_bg*20000
            tmp_count_bg = tmp_count_bg + 1;
            disp(['Current bg event number is : ',num2str(count_gm_bg),' out of ',num2str(tmp_global_count_bg),...
                ' for GMPE ', num2str(GMPEindex(gg))]);
        end
        
          %Perform GMPE for background events 
          %Parameters are GMPE, coefficients, magntiude of event, distance
          %to event, vs30 at site, focal mechanism and depth, fault dip,
          %depths to VS=1 km/s and 2.5 km/s, and vibration periods 
          
          %Parameters are run through GMPE function for each event
          %Returns median ground motion in terms of g with standard deviation for each event 
          
          if GMPEindex(gg)==11 || GMPEindex(gg)==12
          
        repi_bg(:,ii)  = deg2km(distance(site(:,1),site(:,2),bg_catalog_divide(ii,3), bg_catalog_divide(ii,2))); %Determines source to site distance for each event for each site
        rhypo_bg(:,ii)  = sqrt(repi_bg(:,ii) .^2 + fd(ii).^2);          
              
              
           %For vs30 = 300 m/s or USGS vs30 value
        [GMmedtmp,Sintratmp,Sintertmp] = GMPE_Malawi(GMPEindex(gg),COEF,bg_catalog_divide(ii,4),repi_bg(:,ii),rhypo_bg(:,ii),0,repi_bg(:,ii),rhypo_bg(:,ii),vs30_site(:,1),bg_fm(ii,1),fd(ii),Z10_Z25,bg_fm(ii,2),T); 
          %For ref vs30 value
        [GMmedtmp_ref,Sintratmp_ref,Sintertmp_ref] = GMPE_Malawi(GMPEindex(gg),COEF,bg_catalog_divide(ii,4),repi_bg(:,ii),rhypo_bg(:,ii),0,repi_bg(:,ii),rhypo_bg(:,ii),vs30_site(:,2),bg_fm(ii,1),fd(ii),Z10_Z25,bg_fm(ii,2),T); % 
        
        
          else %other GMPE function cannot incorporate multiple sites and multiple SA
              %So run through for loop instead
              
        GMmedtmp=zeros(num_T,num_site); GMmedtmp_ref=zeros(num_T,num_site); 
        Sintratmp=zeros(num_T,num_site); Sintratmp_ref = zeros(num_T,num_site); 
        Sintertmp=zeros(num_T,num_site); Sintertmp_ref= zeros(num_T,num_site);
              
              for ss=1:num_site
              
                repi_bg(ss,ii)  = deg2km(distance(site(ss,1),site(ss,2),bg_catalog_divide(ii,3),bg_catalog_divide(ii,2))); %distance between event and site (map view)
                rhypo_bg(ss,ii) = sqrt(repi_bg(ss,ii).^2 + fd(ii).^2);  %distance between event and site (3D)        
                  
                %For vs30 = 300 m/s or USGS vs30 value
                [GMmedtmp(:,ss),Sintratmp(:,ss),Sintertmp(:,ss)] = GMPE_Malawi(GMPEindex(gg),COEF,bg_catalog_divide(ii,4),repi_bg(ss,ii),rhypo_bg(ss,ii),0,repi_bg(ss,ii),rhypo_bg(ss,ii),vs30_site(ss,1)*ones(num_T,1),bg_fm(ii,1),fd(ii),Z10_Z25,bg_fm(ii,2),T); 
                %For ref vs30 value
                [GMmedtmp_ref(:,ss),Sintratmp_ref(:,ss),Sintertmp_ref(:,ss)] = GMPE_Malawi(GMPEindex(gg),COEF,bg_catalog_divide(ii,4),repi_bg(ss,ii),rhypo_bg(ss,ii),0,repi_bg(ss,ii),rhypo_bg(ss,ii),vs30_site(ss,2)*ones(num_T,1),bg_fm(ii,1),fd(ii),Z10_Z25,bg_fm(ii,2),T); % 
              end
              
        GMmedtmp=GMmedtmp';    GMmedtmp_ref=GMmedtmp_ref'; 
        Sintratmp=Sintratmp'; Sintratmp_ref = Sintratmp_ref'; 
        Sintertmp=Sintertmp'; Sintertmp_ref= Sintertmp_ref';   
        
          end %end if statement
          
        EPS_bg(ii,:,:) = randn(1,num_T,num_site);%Choose normally distributed random number for each SA and site
        
        %Median ground motion for each temporal acceleration multipled by
        %exponent randomly assessed through its standard deviation
        
        for ss=1:num_site
        %Median ground motion for each temporal acceleration multipled by
        %exponent randomly assessed through its standard deviation
      
        GM_bg(gg,ii,:,ss)  = GMmedtmp(ss,:).*exp(sqrt(Sintratmp(ss,:).^2+Sintertmp(ss,:).^2).*(EPS_bg(ii,:,ss)));
        GM_bg_ref(gg,ii,:,ss)  = GMmedtmp_ref(ss,:).*exp(sqrt(Sintratmp_ref(ss,:).^2+Sintertmp_ref(ss,:).^2).*(EPS_bg(ii,:,ss)));
        
        %Take mean GM_bg event for event ii for site ss from all GMPE
            if gg==num_GMPEindex
                mean_GM_bg(ii,:,ss)=squeeze(mean(GM_bg(:,ii,:,ss)))';
                mean_GM_bg_ref(ii,:,ss)=squeeze(mean(GM_bg_ref(:,ii,:,ss)))';
            end    
        end
        
        
    end %end ii loop
      
    clear COEF

end

%Select and store only num_need GM for each SA
for gg = 1:num_GMPEindex
    
    for ss=1:num_site
        
        for ii = 1:num_T
    
        %Sorts background record into ascending order of GM_bg
        tmp = sortrows(GM_bg(gg,:,ii,ss)'); tmp_ref = sortrows(GM_bg_ref(gg,:,ii,ss)');
    
        %Find the num_need_bg number of events
        %For vs30 = USGS value and vs30=ref
        AMAX_GM_bg{gg}{ss}(1:num_need,ii)  = single(tmp(tmp_global_count_bg-num_need+1:tmp_global_count_bg,1)); 
        AMAX_GM_bg_ref{gg}{ss}(1:num_need,ii)  = single(tmp_ref(tmp_global_count_bg-num_need+1:tmp_global_count_bg,1));

        %Store mean_GM for each event, and add catalog data for event (needed for disagg plots)
        if gg==num_GMPEindex
            
            tmp = sortrows([mean_GM_bg(:,ii,ss),EPS_bg(:,ii,ss),rhypo_bg(ss,:)',bg_catalog(divide_bg_indx(1):1:divide_bg_indx(2),4)],1);
            tmp_ref = sortrows([mean_GM_bg_ref(:,ii,ss),EPS_bg(:,ii,ss),rhypo_bg(ss,:)',bg_catalog(divide_bg_indx(1):1:divide_bg_indx(2),4)],1);  
            mean_GM_take{ss}{ii} = single(tmp(tmp_global_count_bg-num_need+1:tmp_global_count_bg,:)); mean_GM_take_ref{ss}{ii}= single(tmp_ref(tmp_global_count_bg-num_need+1:tmp_global_count_bg,:));
        end
        
        end
    end  
end

%% Save parameter based on num_para
%Note stored as single variables to save memory

assignin('base',strcat('AMAX_GM_bg_',string(num_par)),AMAX_GM_bg)
assignin('base',strcat('AMAX_GM_bg_ref_',string(num_par)),AMAX_GM_bg_ref)
assignin('base',strcat('AMAX_GM_mean_bg_',string(num_par)),mean_GM_take)
assignin('base',strcat('AMAX_GM_mean_bg_ref_',string(num_par)),mean_GM_take_ref)

save(strcat('AMAX_GM_bg_',string(num_par)),strcat('AMAX_GM_bg_',string(num_par)),strcat('AMAX_GM_bg_ref_',string(num_par)),...
    strcat('AMAX_GM_mean_bg_',string(num_par)),strcat('AMAX_GM_mean_bg_ref_',string(num_par)),'vs30_site','site','num_site','-v7.3');

%% Selection and creation of combined catalog from Direct and Adapted MSSM 

%load fault geometrical models for sourc2site calc
load fault_geom_1; load fault_geom_2; load sec_geom_1;

%load fault catalogs
load(char(syncat_name),'EQCAT');
load('MSSM_Catalog_Adapted_em','StochasticEventCatalog');

%from all MSSM direct (EQCAT) and/or MSSM adapted take
 
 % 1) event number
 % 2) occurrence time
 % 3) magnitude
 % 4) source_id
 % 5) width case (MSSM adapted only)
%and combine into one cell array

syncat{1}=EQCAT(:,[1,2,5,7]);%MSSM direct catalog

%REMOVE OR COMMENT OUT BEFORE RUNNING ON BC
%syncat{1}=syncat{1}(1:5000,:);
global_count=length(syncat{1});

divide_indx = [floor(global_count*((par_opts(num_par)-1)/length(par_opts)))+1 floor(global_count*(par_opts(num_par)/length(par_opts)))];
syncat{1}=syncat{1}([divide_indx(1):divide_indx(2)],:);

for i=1:length(StochasticEventCatalog)
    syncat{i+1}=StochasticEventCatalog{i}(:,[1,2,4,6,5]);%update with more event catalogs as they are generatred
    
    %REMOVE OF COMMENT OUT BEFORE RUNNING ON BC
    %syncat{i+1}=syncat{i+1}(1:100000,:);
    
    
    global_count=length(syncat{i+1});
    divide_indx = [floor(global_count*((par_opts(num_par)-1)/length(par_opts)))+1 floor(global_count*(par_opts(num_par)/length(par_opts)))];
    syncat{i+1}=syncat{i+1}([divide_indx(1):divide_indx(2)],:);
    
end

num_syncat=length(syncat);

%% Probabilistic seismic hazard analysis for finite fault sources

disp('                                  ');
disp('- - - PSHA for fault sources - - -');

%REMOVE OR COMMENT OUT BEFORE RUNNING ON BC
%num_need=1000;

for ss=1:length(syncat)%loop through each syncat

global_count = length(syncat{ss});
count = 0; group_count = 1;
num_taken = zeros(num_GMPEindex,1); num_taken_ref = zeros(num_GMPEindex,1);

%GM, EPS and GM_ref need to be in cell array length of num_GMPE and
%num_site
AMAX_GM = repmat({repmat({[]},num_site,1)},num_GMPEindex,1); AMAX_GM_ref=repmat({repmat({[]},num_site,1)},num_GMPEindex,1);
AMAX_GM_mean = repmat({[]},num_site,1); AMAX_GM_mean_ref = repmat({[]},num_site,1); 

GM     = single(zeros(num_GMPEindex,num_need,num_T,num_site));
GM_ref = single(zeros(num_GMPEindex,num_need,num_T,num_site));
GM_mean = single(zeros(num_need,num_T,num_site));
GM_mean_ref = single(zeros(num_need,num_T,num_site));
EPS    = single(zeros(num_need,num_T,num_site));
MAG    = zeros(1,num_need);
DIST   = zeros(num_site,4,num_need);
rup_t = zeros(length(global_count),1);
tmp_count = 1;

clock_old = cputime;

for ijk = 1:global_count %For each event in the synthetic record
    
    if ijk == tmp_count*10000
        tmp_count = tmp_count + 1;
        clock_new = cputime;
        disp(['Current event number is : ',num2str(ijk),' out of ',num2str(global_count), ' for catalog ', num2str(ss),....
            ' & Required time is: ',num2str(clock_new-clock_old)]);
        clock_old = clock_new;
    end
    
    %Find grid points for fault that ruptured in catalog       
   syncat_ijk=syncat{ss}(ijk,:);
    
    %directly obtain grids from known earthquake geometry for direct MSSM
    if ss==1
        [gp]=source2site_v1(syncat_ijk,faultgrid_wgs84_1,secgrid_wgs84_1,num_sec,num_fault,grid_indx1);
    
    %float gp on fault grids depending on magnitude for adapted MSSM
    else
        %treat small events as point sources
        if syncat_ijk(3)<5.4
            [gp]=source2site_v2(syncat_ijk,faultgrid_wgs84_1,faultgrid_wgs84_2,grid_indx1,grid_indx2);
            
        %treat large events as 2D planes  
        else
            [gp]=source2site_v3(syncat_ijk,faultgrid_wgs84_1,faultgrid_wgs84_2,num_fault,num_multi_fault,num_MSSM,grid_indx1,grid_indx2);
        end
    end
    
    %Return four types of distances
    %DIST(:,1) = Rjb, 2D distance
    %DIST(:,2) = Rcd, 3D distance
    %DIST(:,3) = Repi, epicentral distance (random point on top of fault plane)
    %DIST(:,4) = Rhypo, hypocentral distance (random point in fault plane) 
    
     distance_tmp=zeros(num_site,4);
     
     for ii=1:num_site
            dist_fault1 = deg2km(distance(site(ii,1),site(ii,2),gp(:,1),gp(:,2))); % Rjb, 2D distance
            dist_fault2 = sqrt(dist_fault1.^2 + gp(:,3).^2); % Rcd, 3D distance
         
           if size(gp,1)==1 %if event treated as a point source, rjb=repi, rcd=rhypo
                distance_tmp(ii,:)=[dist_fault1 dist_fault2 dist_fault1 dist_fault2];
           else 
                ran_dist   = ceil(size(gp,1)*rand(1,1));
                repi  = dist_fault1(ran_dist); %epicentral distance (random point on top of fault plane)
                rhypo = dist_fault2(ran_dist); %hypocentral distance (random point in fault plane) 
                distance_tmp(ii,:)= [min(dist_fault1); min(dist_fault2); repi; rhypo]; 
           end
     end
    
    count=count+1;
    
    %index source dip from MSSM_sources
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
    
    DIST(:,:,count)=distance_tmp;
    MAG(count) = syncat_ijk(3);
    
    %Run through all GMPE for each event in the catalog
    for gg = 1:num_GMPEindex
    
    GMPEID  = GMPEindex(gg);
    assignin('base','COEF',eval(strcat('COEF',num2str(gg))));   
         
    % Ground motion intensity. Run GMPE function for each event, and
    % considers 4 possible distances
  
    if GMPEindex(gg)==11 || GMPEindex(gg)==12
    % USGS vs30 map/vs30 = 300 m/s
    [GMmed,Sintra,Sinter] = GMPE_Malawi(GMPEID,COEF,MAG(count),DIST(:,1,count),DIST(:,2,count),0,DIST(:,3,count),DIST(:,4,count),vs30_site(:,1),FMtype(1),min(gp(:,3)),Z10_Z25,dip,T);
    %Ref vs30 value
    [GMmed_ref,Sintra_ref,Sinter_ref] = GMPE_Malawi(GMPEID,COEF,MAG(count),DIST(:,1,count),DIST(:,2,count),0,DIST(:,3,count),DIST(:,4,count),vs30_site(:,2),FMtype(1),min(gp(:,3)),Z10_Z25,dip,T);
    
    else
  
    GMmed=zeros(num_T,num_site); GMmed_ref=zeros(num_T,num_site); 
    Sintra=zeros(num_T,num_site); Sintra_ref = zeros(num_T,num_site); 
    Sinter=zeros(num_T,num_site); Sinter_ref= zeros(num_T,num_site); 
    
    for ii=1:num_site
  
        % USGS vs30 map/vs30 = 300 m/s
        [GMmed(:,ii),Sintra(:,ii),Sinter(:,ii)] = GMPE_Malawi(GMPEID,COEF,MAG(count),DIST(ii,1,count),DIST(ii,2,count),0,DIST(ii,3,count),DIST(ii,4,count),vs30_site(ii,1)*ones(num_T,1),FMtype(1),min(gp(:,3)),Z10_Z25,dip,T);
        %Ref vs30 value
        [GMmed_ref(:,ii),Sintra_ref(:,ii),Sinter_ref(:,ii)] = GMPE_Malawi(GMPEID,COEF,MAG(count),DIST(ii,1,count),DIST(ii,2,count),0,DIST(ii,3,count),DIST(ii,4,count),vs30_site(ii,2)*ones(num_T,1),FMtype(1),min(gp(:,3)),Z10_Z25,dip,T);
    end
    
    GMmed=GMmed';   GMmed_ref=GMmed_ref'; 
    Sintra=Sintra'; Sintra_ref = Sintra_ref'; 
    Sinter=Sinter'; Sinter_ref= Sinter_ref';  
    
    end %end if statement
    
     EPS(count,:,:) = randn(1,num_T,num_site);

    for ii=1:num_site
    
        GM(gg,count,:,ii) = GMmed(ii,:).*exp(sqrt(Sintra(ii,:).^2+Sinter(ii,:).^2).*EPS(count,:,ii));% Ground motion for each event randomly altered by standard deviation
        GM_ref(gg,count,:,ii) = GMmed_ref(ii,:).*exp(sqrt(Sintra_ref(ii,:).^2+Sinter_ref(ii,:).^2).*EPS(count,:,ii));% Ground motion for each event randomly altered by standard deviation
    
        %mean GM for each GMPE
    if gg==num_GMPEindex
        GM_mean(count,:,ii)=mean(GM(:,count,:,ii));
        GM_mean_ref(count,:,ii)=mean(GM_ref(:,count,:,ii));
    end
        
    end
    
    if ijk == group_count*num_need || ijk==global_count %Number of events in group size
        
        for ii=1:num_site
        
        if group_count==1
                AMAX_GM{gg}{ii} = squeeze(GM(gg,:,:,ii)); 
                AMAX_GM_ref{gg}{ii} = squeeze(GM_ref(gg,:,:,ii));
                
                if gg==num_GMPEindex
                    for kk=1:num_T
                        AMAX_GM_mean{ii}(1:num_need,1:4,kk)=[GM_mean(:,kk,ii),EPS(:,kk,ii),squeeze(DIST(ii,4,:)),MAG(1,:)'];
                        AMAX_GM_mean_ref{ii}(1:num_need,1:4,kk)=[GM_mean_ref(:,kk,ii),EPS(:,kk,ii),squeeze(DIST(ii,4,:)),MAG(1,:)'];
                    end
                end
        else
               for kk = 1:num_T
                    %combine with AMAX_GM from previous groups
                    tmp = [AMAX_GM{gg}{ii}(:,kk); squeeze(GM(gg,:,kk,ii))'];
                    tmp = sortrows(tmp,1,'descend'); %Sort rows in descending order of GM
                
                    tmp_ref = [AMAX_GM_ref{gg}{ii}(:,kk); squeeze(GM_ref(gg,:,kk,ii))'];
                    tmp_ref = sortrows(tmp_ref,1,'descend'); %Sort rows in descending order of GM
                    %only take highest ground motions from most recent simulations (max_GM)
                    %and previous runs (AMAX_GM)
                    AMAX_GM{gg}{ii}(1:num_need,kk) = tmp(1:num_need);
                    AMAX_GM_ref{gg}{ii}(1:num_need,kk) = tmp_ref(1:num_need);
                    
                    if gg==num_GMPEindex
                        %sample EPS, rhypo, and magnitude for deagg analysis
                        tmp_store=vertcat(AMAX_GM_mean{ii}(:,:,kk),[GM_mean(:,kk,ii),EPS(:,kk,ii),squeeze(DIST(ii,4,:)),MAG(1,:)']);
                        tmp_store = sortrows(tmp_store,1,'descend'); %Sort rows in descending order of GM
                        tmp_store_ref=vertcat(AMAX_GM_mean_ref{ii}(:,:,kk),[GM_mean_ref(:,kk,ii),EPS(:,kk,ii),squeeze(DIST(ii,4,:)),MAG(1,:)']);
                        tmp_store_ref = sortrows(tmp_store_ref,1,'descend');
                        AMAX_GM_mean{ii}(1:num_need,:,kk)= tmp_store(1:num_need,1:4);
                        AMAX_GM_mean_ref{ii}(1:num_need,:,kk)= tmp_store_ref(1:num_need,1:4);
                    end
                end
                
        end %end if statement for group count==1
        end %end for loop for num_site
       
    end %end groupcount if statement  
    
    end %end gg loop for GMPE
    
        if ijk == group_count*num_need 
         
        group_count = group_count + 1;
        count=0; 
         
        %Creates new vector of zero now that group is full
        GM     = single(zeros(num_GMPEindex,num_need,num_T,num_site));
        GM_ref = single(zeros(num_GMPEindex,num_need,num_T,num_site));
        GM_mean = single(zeros(num_need,num_T,num_site));
        GM_mean_ref = single(zeros(num_need,num_T,num_site));
        EPS    = single(zeros(num_need,num_T,num_site));
        MAG    = zeros(1,num_need);
        DIST   = zeros(num_site,4,num_need);
        
        end
    
end %end ijk loop for each event

%Save ground motions where each cell represents source ss
%Within each cell array further cells represents ground motions from GMPE gg
AMAX_GM_em{ss}=AMAX_GM;
AMAX_GM_ref_em{ss}=AMAX_GM_ref;

AMAX_GM_mean_tmp{ss}=AMAX_GM_mean;
AMAX_GM_mean_tmp_ref{ss}=AMAX_GM_mean_ref;

end

AMAX_GM_mean_em=cell(num_site,1);
AMAX_GM_mean_em_ref=cell(num_site,1);

%Combine GM from different syncats
for ii=1:num_site
    for kk=1:num_T
        
         tmp_sort=zeros(1,4); tmp_sort_ref=zeros(1,4);
    
       for ss=1:num_syncat
          tmp_sort=vertcat(tmp_sort,AMAX_GM_mean_tmp{ss}{ii}(:,:,kk)); 
          tmp_sort_ref=vertcat(tmp_sort_ref,AMAX_GM_mean_tmp_ref{ss}{ii}(:,:,kk)); 
       end
        AMAX_GM_mean_em{ii}{kk}=single(tmp_sort);  AMAX_GM_mean_em_ref{ii}{kk}=single(tmp_sort_ref); 
    end
    
end


%% Save parameter based on num_para
%Note stored as single variables to save memory
assignin('base',strcat('AMAX_GM_',string(num_par)),AMAX_GM_em);
assignin('base',strcat('AMAX_GM_ref_',string(num_par)),AMAX_GM_ref_em);
assignin('base',strcat('AMAX_GM_mean_',string(num_par)),AMAX_GM_mean_em);
assignin('base',strcat('AMAX_GM_mean_ref_',string(num_par)),AMAX_GM_mean_em_ref);

run_time=toc;

save(strcat('AMAX_GM_',string(num_par)),strcat('AMAX_GM_',string(num_par)),strcat('AMAX_GM_ref_',string(num_par)),...
    strcat('AMAX_GM_mean_',string(num_par)), strcat('AMAX_GM_mean_ref_',string(num_par)),'run_time','-v7.3');
