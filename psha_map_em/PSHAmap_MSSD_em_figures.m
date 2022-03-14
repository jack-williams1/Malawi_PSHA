%% Combines simulation results and makes PSHA Map

mydir  = pwd; idcs   = strfind(mydir,'/');
addpath(mydir(1:idcs(end)-1));

%add path to GM from Blue Crystal Run Updated as necessary
addpath([mydir(1:idcs(end-1)-1),'/psha_stored_gm/20220302_GM']);


load('syncat_PSHA_MSSD_input','NumSimu','vs30_site_ref','prob_level','obs_duration','num_GMPEindex',...
    'lower_return_period','num_need','syncat_name','par_opts','Region','PSHA_Zone',...
    'site_corner_SW','site_corner_NE','site_grid_interval','T_map','site','site_name');

site_spec=site;%sites for site specific psha
clear site
lat_site_vec = site_corner_SW(1):site_grid_interval:site_corner_NE(1);
lon_site_vec = site_corner_SW(2):site_grid_interval:site_corner_NE(2);

[Lat_site,Lon_site] = meshgrid(lat_site_vec,lon_site_vec);

site = [Lat_site(:) Lon_site(:)];

num_site = size(site,1); num_syncat=5;

%if code run before, load sorted and clean ground motions and skip straight
%to 'Seismic Hazard fractiles'
load GM_MSSD_em_20220302

%% load and combine results from each paralleisation

tmp_GM_bg=cell(num_GMPEindex,1); tmp_GM_bg_ref=cell(num_GMPEindex,1);
tmp_GM_fs=cell(num_syncat,1); tmp_GM_fs_ref=cell(num_syncat,1);

count_site=0; tmp_count_site=1;

GM_bg = single(zeros(num_need,num_site,num_GMPEindex)); 
GM_bg_ref = single(zeros(num_need,num_site,num_GMPEindex)); 
GM_fs = repmat({single(zeros(num_need,num_site,num_GMPEindex))},num_syncat,1); 
GM_fs_ref = repmat({single(zeros(num_need,num_site,num_GMPEindex))},num_syncat,1); 


for pp=1:length(par_opts)
    
    disp(['on parallelization ',num2str(pp) ' out of ', num2str(length(par_opts))])
    
   for jj=1:num_GMPEindex 
       
    if floor(pp/2)~=pp/2 
       tmp_GM_bg{jj}=single(zeros(1,num_site)); tmp_GM_bg_ref{jj}=single(zeros(1,num_site));
    end
    
   load(strcat('AMAX_GM_bg_',string(par_opts(pp))))
   
   tmp_GM_bg{jj}=vertcat(tmp_GM_bg{jj}, eval(strcat('AMAX_GM_bg_em_',string(par_opts(pp)),'{',string(num2str(jj)),'}')));
   tmp_GM_bg_ref{jj}=vertcat(tmp_GM_bg_ref{jj}, eval(strcat('AMAX_GM_bg_em_ref_',string(par_opts(pp)),'{',string(num2str(jj)),'}')));
   
   end
   
   load(strcat('AMAX_GM_fs_',string(par_opts(pp))))
   
   for ss=1:num_syncat %for each MSSD syncat
       if floor(pp/2)~=pp/2 
           tmp_GM_fs{ss}={}; tmp_GM_fs_ref{ss}={}; 
           tmp_GM_fs{ss}=single(zeros(1,num_site,num_GMPEindex)); tmp_GM_fs_ref{ss}=single(zeros(1,num_site,num_GMPEindex));
       end
    tmp_GM_fs{ss}=vertcat(tmp_GM_fs{ss}, eval(strcat('AMAX_GM_fs_',string(par_opts(pp)),'{',string(num2str(ss)),'}')));
    tmp_GM_fs_ref{ss}=vertcat(tmp_GM_fs_ref{ss}, eval(strcat('AMAX_GM_fs_ref_',string(par_opts(pp)),'{',string(num2str(ss)),'}')));
   end
   
   clear(strcat('AMAX_GM_bg_em_',string(par_opts(pp)))); clear(strcat('AMAX_GM_bg_em_ref_',string(par_opts(pp))));
   clear(strcat('AMAX_GM_fs_',string(par_opts(pp)))); clear(strcat('AMAX_GM_fs_ref_',string(par_opts(pp))));
   
   %to stop files getting to big sort, find num_need samples for every
   %other loop
   
   if floor(pp/2)==pp/2
    
        for jj=1:num_GMPEindex
           tmp=sort(vertcat(tmp_GM_bg{jj},GM_bg(:,:,jj)),'ascend'); 
           GM_bg(:,:,jj) = single(tmp(length(tmp)-num_need+1:length(tmp),:));

           tmp=sort(vertcat(tmp_GM_bg_ref{jj},GM_bg_ref(:,:,jj)),'ascend');
           GM_bg_ref(:,:,jj) = single(tmp(length(tmp)-num_need+1:length(tmp),:));
        
            for ss=1:num_syncat
        
                tmp=sort(vertcat(tmp_GM_fs{ss}(:,:,jj),GM_fs{ss}(:,:,jj)),'ascend');
                GM_fs{ss}(:,:,jj) = single(tmp(length(tmp)-num_need+1:length(tmp),:));
                
                tmp=sort(vertcat(tmp_GM_fs_ref{ss}(:,:,jj),GM_fs_ref{ss}(:,:,jj)),'ascend');
                GM_fs_ref{ss}(:,:,jj) = single(tmp(length(tmp)-num_need+1:length(tmp),:));

            end
            
        end %end jj loop
   
   tmp_GM_bg=cell(num_GMPEindex,1); tmp_GM_bg_ref=cell(num_GMPEindex,1);
   tmp_GM_fs=cell(num_syncat,1); tmp_GM_fs_ref=cell(num_syncat,1);
   
   end %end if statement
 
end

%% Combine bg and fs ground motions

GM_cb = repmat({single(zeros(num_need,num_site,num_GMPEindex))},num_syncat,1); 
GM_cb_ref = repmat({single(zeros(num_need,num_site,num_GMPEindex))},num_syncat,1); 

for ss=1:num_syncat
    for jj=1:num_GMPEindex
                %combine fault and areal sources and select only num_need
                tmp_cb=sort(vertcat(GM_fs{ss}(:,:,jj),GM_bg(:,:,jj))); tmp_cb=sort(tmp_cb,'ascend');
                GM_cb{ss}(:,:,jj)=tmp_cb(length(tmp_cb)-num_need+1:length(tmp_cb),:);
        
                tmp_cb=sort(vertcat(GM_fs_ref{ss}(:,:,jj),GM_bg_ref(:,:,jj))); tmp_cb=sort(tmp_cb,'ascend');
                GM_cb_ref{ss}(:,:,jj)=tmp_cb(length(tmp_cb)-num_need+1:length(tmp_cb),:);
    end
end

%save('GM_MSSD_em_20220302','GM_bg','GM_bg_ref','GM_fs','GM_fs_ref','GM_cb','GM_cb_ref','-v7.3');

%% Seismic hazard map fractiles

PSA_fractile_CDF_bg = zeros(length(prob_level),num_site,num_GMPEindex); PSA_fractile_CDF_bg_ref = zeros(length(prob_level),num_site,num_GMPEindex);
PSA_fractile_CDF_fs = zeros(length(prob_level),num_site,num_syncat*num_GMPEindex); PSA_fractile_CDF_fs_ref = zeros(length(prob_level),num_site,num_syncat*num_GMPEindex);
PSA_fractile_CDF_cb = zeros(length(prob_level),num_site,num_syncat*num_GMPEindex); PSA_fractile_CDF_cb_ref = zeros(length(prob_level),num_site,num_syncat*num_GMPEindex);

PSA_fractile_mean_bg=zeros(length(prob_level),num_site); PSA_fractile_mean_bg_ref=zeros(length(prob_level),num_site);
PSA_fractile_mean_fs=zeros(length(prob_level),num_site); PSA_fractile_mean_fs_ref=zeros(length(prob_level),num_site);
PSA_fractile_mean=zeros(length(prob_level),num_site); PSA_fractile_cov=zeros(length(prob_level),num_site);
PSA_fractile_mean_ref=zeros(length(prob_level),num_site); PSA_fractile_cov_ref=zeros(length(prob_level),num_site);
PSA_fractile_iqr=zeros(length(prob_level),num_site); PSA_fractile_iqr_ref=zeros(length(prob_level),num_site); 

psha_comp=zeros(length(prob_level)*2,2,num_site); psha_comp_ref=zeros(length(prob_level)*2,2,num_site);


for ii = 1:length(prob_level)
    
    count=0;
    for ss =1:num_syncat
        
        %Return 20 realisations of hazard levels from all GMM-source model combinations
        for jj=1:num_GMPEindex
            count=count+1;
            PSA_fractile_CDF_fs(ii,:,count) = GM_fs{ss}(num_need-round(obs_duration/prob_level(ii))+1,:,jj);
            PSA_fractile_CDF_fs_ref(ii,:,count) = GM_fs_ref{ss}(num_need-round(obs_duration/prob_level(ii))+1,:,jj);
            PSA_fractile_CDF_cb(ii,:,count) = GM_cb{ss}(num_need-round(obs_duration/prob_level(ii))+1,:,jj);
            PSA_fractile_CDF_cb_ref(ii,:,count) = GM_cb_ref{ss}(num_need-round(obs_duration/prob_level(ii))+1,:,jj);
            
            if ss==1
                PSA_fractile_CDF_bg(ii,:,jj) = GM_bg(num_need-round(obs_duration/prob_level(ii))+1,:,jj);
                PSA_fractile_CDF_bg_ref(ii,:,jj) = GM_bg_ref(num_need-round(obs_duration/prob_level(ii))+1,:,jj);
            end
        end
    end
    
    %Difference between G-R and (1) Direct MSSD and (2) Char ground motions
    %For char and G-R length limited cases, and Boore 2014 GMM   
    psha_comp(ii,:,:)=vertcat(PSA_fractile_CDF_cb(ii,:,9)-PSA_fractile_CDF_cb(ii,:,1),...
    PSA_fractile_CDF_cb(ii,:,9)-PSA_fractile_CDF_cb(ii,:,5));

    psha_comp_ref(ii,:,:)=vertcat(PSA_fractile_CDF_cb_ref(ii,:,9)-PSA_fractile_CDF_cb_ref(ii,:,1),...
    PSA_fractile_CDF_cb_ref(ii,:,9)-PSA_fractile_CDF_cb_ref(ii,:,5));

    %Difference between layer and length limited widths for G-R model and Boore 2014 GMM
    psha_comp(ii+2,1,:)=PSA_fractile_CDF_cb(ii,:,17)-PSA_fractile_CDF_cb(ii,:,9);
    psha_comp_ref(ii+2,1,:)=PSA_fractile_CDF_cb_ref(ii,:,17)-PSA_fractile_CDF_cb_ref(ii,:,9);
    
    %Derive mean and cov values for each site
    for tt=1:num_site
        PSA_fractile_mean_bg(ii,tt)=mean(PSA_fractile_CDF_bg(ii,tt,:));
        PSA_fractile_mean_fs(ii,tt)=mean(PSA_fractile_CDF_fs(ii,tt,:));
        PSA_fractile_mean(ii,tt)=mean(PSA_fractile_CDF_cb(ii,tt,:));
        std_tmp=std(PSA_fractile_CDF_cb(ii,tt,:));
        PSA_fractile_cov(ii,tt)=std_tmp/PSA_fractile_mean(ii,tt);
        PSA_fractile_iqr(ii,tt)=iqr(PSA_fractile_CDF_cb(ii,tt,:));
        
        PSA_fractile_mean_bg_ref(ii,tt)=mean(PSA_fractile_CDF_bg_ref(ii,tt,:));
        PSA_fractile_mean_fs_ref(ii,tt)=mean(PSA_fractile_CDF_fs_ref(ii,tt,:));
        PSA_fractile_mean_ref(ii,tt)=mean(PSA_fractile_CDF_cb_ref(ii,tt,:));
        std_tmp=std(PSA_fractile_CDF_cb_ref(ii,tt,:));
        PSA_fractile_cov_ref(ii,tt)=std_tmp/PSA_fractile_mean_ref(ii,tt);
        PSA_fractile_iqr_ref(ii,tt)=iqr(PSA_fractile_CDF_cb_ref(ii,tt,:));
    end
    
end

    
%% Load data for plotting maps

lat_site_vec = site_corner_SW(1):site_grid_interval:site_corner_NE(1);
lon_site_vec = site_corner_SW(2):site_grid_interval:site_corner_NE(2);

[Lat_site,Lon_site] = meshgrid(lat_site_vec,lon_site_vec);

site = [Lat_site(:) Lon_site(:)];

addpath([mydir(1:idcs(end)-1) '/gis_files']);
addpath([mydir(1:idcs(end)-1) '/misc_functions']);

ellipsoid = wgs84Ellipsoid('meters');
utms = defaultm('utm'); utms.zone = '36L'; utms.geoid = ellipsoid; 
utms.flatlimit = []; utms.maplatlimit = []; utms = defaultm(utms);

load map_data_EastAfrica
load LakeMalawiBorder

LakeMalawi = shaperead('malawi_lake.shp');
LakeMalawiCoord = [LakeMalawi.Y(1,1:end-1)' LakeMalawi.X(1,1:end-1)'];

MSSD = shaperead('MSSD_fault.shp');%Update if MSSD is revised
num_traces = length(MSSD);

cmap = crameri('batlow');

%% Surface plot of ground shaking level for each site

%vs30 value to use for PSHA maps
%1 equals USGS map, else is ref vs30 value
opt = 2;

if opt==1
    PSA_fractile_mean_bg_opt = PSA_fractile_mean_bg;
    PSA_fractile_mean_fs_opt = PSA_fractile_mean_fs;
    PSA_fractile_mean_opt = PSA_fractile_mean;
    PSA_fractile_cov_opt = PSA_fractile_cov;
    PSA_fractile_iqr_opt = PSA_fractile_iqr;
    psha_comp_opt=psha_comp; 
else
    PSA_fractile_mean_bg_opt = PSA_fractile_mean_bg_ref;
    PSA_fractile_mean_fs_opt = PSA_fractile_mean_fs_ref;
    PSA_fractile_mean_opt = PSA_fractile_mean_ref;
    PSA_fractile_cov_opt = PSA_fractile_cov_ref;
    PSA_fractile_iqr_opt = PSA_fractile_iqr_ref;
    psha_comp_opt=psha_comp_ref; 
end

clear tmp1_ tmp2_ tmp3_



%% Plot seismic hazard maps

figure (1000)
tiledlayout(2,3,'tilespacing','compact')

plabel_opt=strcat(["10% PoE in 50 years","2% PoE in 50 years"]);
label_opt=vertcat(strcat(["(a) Areal Sources","(b) MSSD Sources","(c) Combined"]),strcat(["(d) Areal Sources","(e) MSSD Sources","(f) Combined"]));


for ii = 1:length(prob_level)
    
    tmp1_{ii} = reshape(PSA_fractile_mean_bg_opt(ii,:),[length(lon_site_vec),length(lat_site_vec)]);
    tmp2_{ii} = reshape(PSA_fractile_mean_fs_opt(ii,:),[length(lon_site_vec),length(lat_site_vec)]);
    tmp3_{ii} = reshape(PSA_fractile_mean_opt(ii,:),[length(lon_site_vec),length(lat_site_vec)]);
     
    nexttile
    surf(Lon_site,Lat_site,tmp1_{ii}); hold on; shading interp; view(0,90); axis equal; axis(Region);     
    title(label_opt(ii,1),'Fontweight','normal'); subtitle([plabel_opt(ii),'']) ;
    hold on; caxis([0.0 0.8]); set(gca,'fontsize',9,'TitleFontSizeMultiplier',1.365);
    fill3(LakeMalawiCoord(:,2),LakeMalawiCoord(:,1),1000*ones(length(LakeMalawiCoord),1),'w'); hold on
    plot3(MapData2(:,2),MapData2(:,1),1000*ones(length(MapData2(:,1)),1),'w-');
    plot3(LakeMalawiBorder(:,2),LakeMalawiBorder(:,1),1000*ones(length(LakeMalawiBorder),1),'k-');hold on
    plot3(site_spec(:,2),site_spec(:,1),1000*ones(length(site_spec),1),'ks','markerfacecolor','w','markersize',7);
    
    for ss=1:length(site_spec)
        text(site_spec(ss,2)-0.8,site_spec(ss,1)-0.2,1000,site_name(ss),'Color', 'w','FontSize',8.5);
    end
    
    
    %subplot(length(prob_level),3,(3*(ii-1))+2);
    nexttile
    surf(Lon_site,Lat_site,tmp2_{ii}); hold on; shading interp; view(0,90); axis equal; axis(Region); title(label_opt(ii,2),'Fontweight','normal');         
    subtitle([plabel_opt(ii),'']); 
    fill3(LakeMalawiCoord(:,2),LakeMalawiCoord(:,1),1000*ones(length(LakeMalawiCoord),1),'w');
    plot3(MapData2(:,2),MapData2(:,1),1000*ones(length(MapData2(:,1)),1),'w-');
    plot3(LakeMalawiBorder(:,2),LakeMalawiBorder(:,1),1000*ones(length(LakeMalawiBorder),1),'k-');hold on
    hold on; caxis([0.0 0.8]); set(gca,'XTick',[], 'YTick', []); set(gca,'fontsize',11);
    for jj = 1:num_traces    
        [LAT,LON] = minvtran(utms,MSSD(jj).X,MSSD(jj).Y);
        plot3(LON,LAT,1000*ones(length(LON),1),'r','linewidth',0.75); hold on;
    end
    
    plot3(site_spec(:,2),site_spec(:,1),1000*ones(length(site_spec),1),'ks','markerfacecolor','w','markersize',7);

    %subplot(length(prob_level),3,(3*(ii-1))+3);
    nexttile
    surf(Lon_site,Lat_site,tmp3_{ii}); hold on; shading interp; view(0,90); axis equal; axis(Region); title(label_opt(ii,3),'Fontweight','normal');        
    subtitle([plabel_opt(ii),'']);
    fill3(LakeMalawiCoord(:,2),LakeMalawiCoord(:,1),1000*ones(length(LakeMalawiCoord),1),'w');
    plot3(MapData2(:,2),MapData2(:,1),1000*ones(length(MapData2(:,1)),1),'w-');
    plot3(LakeMalawiBorder(:,2),LakeMalawiBorder(:,1),1000*ones(length(LakeMalawiBorder),1),'k-');hold on
    set(gca,'XTick',[], 'YTick', [])
    hold on; caxis([0.0 0.8]); colormap(cmap); h5=colorbar; h5.Label.String='PGA (g)'; set(gca,'fontsize',11);
    for jj = 1:num_traces    
        [LAT,LON] = minvtran(utms,MSSD(jj).X,MSSD(jj).Y);
        plot3(LON,LAT,1000*ones(length(LON),1),'r','linewidth',0.75); hold on;
         
    end
   plot3(site_spec(:,2),site_spec(:,1),1000*ones(length(site_spec),1),'ks','markerfacecolor','w','markersize',7);
    
   %clear tmp1 tmp2 tmp3
    
end

set(gcf,'position',[395 91 485 706])


%% Plot figure for seismic hazard spatial distribution differences uncertainity

figure (1001)

%select PoE level to run figure for by uncommenting
plevel=1;% 10% PoE in 50 years
%plevel=2;% 2% PoE in 50 years

label_opt=vertcat(strcat(["(a) 10% PoE in 50 years", "(b) Interquartile Range","(c) CoV"]),strcat(["(d) GR-MSSD Direct", "(e) GR-Char","(f) Layer Limited-Length Limited"]));

tiledlayout(2,3,'tilespacing','tight')

for ii = 1:3
    
    if ii==1
    tmp6 = tmp3_{plevel};   
    tmp6a=([0 0.4]); tmp6b='PGA (g)'; tmp6c='r';
    
    tmp7 = reshape(psha_comp_opt(plevel,ii,:),[length(lon_site_vec),length(lat_site_vec)]);
    
    elseif ii ==2
    tmp6 = reshape(PSA_fractile_iqr_opt(plevel,:),[length(lon_site_vec),length(lat_site_vec)]);
    tmp6a=([0 0.2]); tmp6b='Interquartile range (g)'; tmp6c='k';   
        
    tmp7 = reshape(psha_comp_opt(plevel,ii,:),[length(lon_site_vec),length(lat_site_vec)]);
    
    elseif ii==3
    tmp6 = reshape(PSA_fractile_cov_opt(plevel,:),[length(lon_site_vec),length(lat_site_vec)]);
    tmp6a=([0 1]); tmp6b= 'CoV' ;tmp6c='k';
    
    tmp7 = reshape(psha_comp_opt(ii,1,:),[length(lon_site_vec),length(lat_site_vec)]);
    
    end
    
    %Plot (a), (c) and (e)
    nexttile(ii)
    if ii==1
        cmap1 = crameri('batlow');
    else
        cmap1 = 'pink';
    end
    
    surf(Lon_site,Lat_site,tmp6); hold on; shading interp; view(0,90); axis equal; axis(Region);     
    title(label_opt(1,ii),'Fontweight','normal'); hold on; caxis(tmp6a); 
    fill3(LakeMalawiCoord(:,2),LakeMalawiCoord(:,1),1000*ones(length(LakeMalawiCoord),1),'w'); hold on
    plot3(MapData2(:,2),MapData2(:,1),1000*ones(length(MapData2(:,1)),1),'w-');
    plot3(LakeMalawiBorder(:,2),LakeMalawiBorder(:,1),1000*ones(length(LakeMalawiBorder),1),'k-');hold on
    colormap(gca,cmap1); h4=colorbar; h4.Label.String=tmp6b;
    
    for jj = 1:num_traces    
        [LAT,LON] = minvtran(utms,MSSD(jj).X,MSSD(jj).Y);
        plot3(LON,LAT,1000*ones(length(LON),1),tmp6c,'linewidth',0.75); hold on;
    end
    
    if ii==1
        set(gca,'fontsize',11);
        plot3(site_spec(:,2),site_spec(:,1),1000*ones(length(site_spec),1),'ks','markerfacecolor','w','markersize',7);
        tmp8=vertcat([0.8 -0.15 0.55],[0.2 0.2 0.2]);
        for ss=1:length(site_spec)
            text(site_spec(ss,2)-tmp8(1,ss),site_spec(ss,1)+tmp8(2,ss),1000,site_name(ss),'Color', 'w','FontSize',9);
        end
    else
        set(gca,'XTick',[], 'YTick', [],'fontsize',11,'TitleFontSizeMultiplier',1.09);
    end
   
    box on;
    hold off
    
    %Plot (b), (d) and (f)
    nexttile(ii+3)
    cmap = crameri('vik');
    surf(Lon_site,Lat_site,tmp7); hold on; shading interp; view(0,90); axis equal; axis(Region); title(label_opt(2,ii),'Fontweight','normal');         
    fill3(LakeMalawiCoord(:,2),LakeMalawiCoord(:,1),1000*ones(length(LakeMalawiCoord),1),'w');
    plot3(MapData2(:,2),MapData2(:,1),1000*ones(length(MapData2(:,1)),1),'k-');
    hold on; caxis([-0.4 0.4]); colormap(gca,cmap);  set(gca,'XTick',[], 'YTick', [],'fontsize',11,'TitleFontSizeMultiplier',1.09);
    for jj = 1:num_traces    
        [LAT,LON] = minvtran(utms,MSSD(jj).X,MSSD(jj).Y);
        plot3(LON,LAT,1000*ones(length(LON),1),'r','linewidth',0.75); hold on;
    end

    if ii==3
        h4=colorbar; h4.Label.String='PGA Difference (g)';
    end
      box on;
    hold off


end

set(gcf,'position',[456 132 777 845]);


%% Write 10% PoE values to csv file

lat = string(lat_site_vec);
lon = string(lon_site_vec);

prob_level=[475 2475];

for ii=1:length(prob_level)

pga=tmp3_{ii}';%combined PGA values for given prob_level

MalawiPGA = array2table(pga,'RowNames',lat,'VariableNames',lon);
MalawiPGA = flipud(MalawiPGA); % So northern most latitudes appear at top of table

if opt==1
writetable(MalawiPGA,strcat('MalawiPGA_USGS_vs30_',string(prob_level(ii))),'WriteRowNames',true);
else
writetable(MalawiPGA,strcat('MalawiPGA_ref_vs30_',string(prob_level(ii))),'WriteRowNames',true);
end

end