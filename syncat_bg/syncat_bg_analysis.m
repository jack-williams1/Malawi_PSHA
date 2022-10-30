%% Syncat bg analysis

load syncat_bg
load new_SZ
malawi_bg_sz = [new_SZ(7) new_SZ(8) new_SZ(9) new_SZ(14) new_SZ(16) new_SZ(17) new_SZ(19) new_SZ(20) new_SZ(23) new_SZ(26) new_SZ(27)];

mydir  = pwd; idcs   = strfind(mydir,'/');
addpath(mydir(1:idcs(end)-1));
addpath([mydir(1:idcs(end)-1) '/gis_files']);

load('syncat_PSHA_MSSM_input','NumSimu','bg_mmin','PSHA_Zone');

%% QC: Find moment rate in bg_catalog for events in source zone

sz =[2 3 5 9 10 11]; %All source zones within Malawi to test based on index in bg_source 
%set to 3 to test for Malawi

for ss=1:length(sz)

index_sz=find(bg_catalog(:,6)==sz(ss) & bg_catalog(:,1)<0);%Find events from that source zone and in source area

mo_sz=zeros(length(index_sz),1);

    for k=1:length(index_sz)
        mo_sz(k)= 10.^((bg_catalog(index_sz(k),4)*1.5)+9.05);
    end

mo_sz_r(ss)=sum(mo_sz)/NumSimu;

end

%% Derive mag-freq relation for source zone from discretized G-R curve
%Uses a value scaled for events in Malawi from seismic_moment_spreadsheet
%Assumes truncated exponential relation

bg_a = [2.9 4.7 2.8 2.6 0.4 0.1]; %see inMalawiG-R spreadsheet
GR_Discrete={[] [] [] [] [] []};  syncat_bg_AnnualRate={[] [] [] [] [] []};

for ss=1:length(sz)

    bg_b=malawi_bg_sz(sz(ss)).b;
    mmin = bg_mmin;
    dm = 0.01;

    %testing sz3 for higher Mwmax
    %{
    if ss==2
    malawi_bg_sz(sz(ss)).Mwmax=8.1;
    end
    %}
    
    bg_mag_range_GR = (4.5:dm:malawi_bg_sz(sz(ss)).Mwmax);

    bg_AnnualRate=zeros(length(bg_mag_range_GR),1);

    beta = log(10)*malawi_bg_sz(sz(ss)).b;
    rate = 10^(bg_a(ss)-bg_b*mmin);
    GR_Discrete{ss}=zeros(length(bg_mag_range_GR),6);
    
    syncat_bg_AnnualRate{ss}=zeros(length(bg_mag_range_GR),1);
    index_sz=find(bg_catalog(:,6)==sz(ss) & bg_catalog(:,1)<0);
    
    for ii = 1:length(bg_mag_range_GR)
        
        GR_Discrete{ss}(ii,1) = bg_mag_range_GR(ii);
        %cum distribution function for non exceedance of m>M for each mag bin
        %(i.e function = 1 for Mmax, 0 for mmin)
        GR_Discrete{ss}(ii,2) = (1 - exp(-beta*(bg_mag_range_GR(ii)-(dm/2) - mmin)))/(1 - exp(-beta*(malawi_bg_sz(sz(ss)).Mwmax-mmin)));
        GR_Discrete{ss}(ii,3) = (1 - exp(-beta*(bg_mag_range_GR(ii)+(dm/2) - mmin)))/(1 - exp(-beta*(malawi_bg_sz(sz(ss)).Mwmax-mmin)));
  
        %probability for m = M
        GR_Discrete{ss}(ii,4) = GR_Discrete{ss}(ii,3) - GR_Discrete{ss}(ii,2);
  
        %cum distribution function for Exceedance of m>M for each mag bin
        %(i.e function = 0 for Mmax, 1 for mmin), multipled by rate
        GR_Discrete{ss}(ii,5) = rate*(1-(1 - exp(-beta*(bg_mag_range_GR(ii) - mmin)))/(1 - exp(-beta*(malawi_bg_sz(sz(ss)).Mwmax-mmin))));
  
         %Moment rate for each mag bin
        GR_Discrete{ss}(ii,6) = GR_Discrete{ss}(ii,4)*rate*10.^((bg_mag_range_GR(ii)*1.5)+9.05);
        
        
        syncat_bg_AnnualRate{ss}(ii) = length(find(bg_catalog(index_sz,4) >= bg_mag_range_GR(ii)))/(NumSimu);    
   end

mo_sz_r1(ss)=sum(GR_Discrete{ss}(:,6)); %moment rate from discretising G-R
end

%% Plot mag freq distribution
figure(23);

for ss=1:length(sz)

semilogy(GR_Discrete{ss}(:,1),GR_Discrete{ss}(:,5),'k-','LineWidth',1.5); hold on; %MFD from theoritical G-R curve
semilogy(GR_Discrete{ss}(:,1),syncat_bg_AnnualRate{ss}(:),'r-','LineWidth',1.5); hold on; %MFD from syncat

end

set(gca,'fontsize',13);  axis([5 8 10^-5 10^0]); legend('Discretized G-R','Event Catalog','Location','northeast');hold on;
xlabel('Magnitude'); ylabel('Annual frequency of exceedance'); grid on; axis square;

%% Plot event locations (Optional)
%plots location of events from a 1000 year record of events
load syncat_bg_tmp

Region2 = [29.6 39 -20.2 -6]; %For plots at larger scale
load map_data_EastAfrica
buffer_200=shaperead('PSHA_BufferArea200_wgs.shp');

index_loc={};

for i=1:length(malawi_bg_sz)
    index_loc{i}=find(bg_catalog(:,6)==i,100);
end

figure(102);
%Need to check figure is plotting correctly for source zones with holes in
fill(new_SZ(7).ZONE(:,1),new_SZ(7).ZONE(:,2),[0.7 0.6 0.2],'FaceAlpha',0.2); hold on; %Rovuma Basin-1
fill(new_SZ(8).ZONE(:,1),new_SZ(8).ZONE(:,2),'r','FaceAlpha',0.2); hold on; %Western Rift-Tanganyika
fill(new_SZ(9).ZONE(:,1),new_SZ(9).ZONE(:,2),'b','FaceAlpha',0.2); hold on; % Malawi Rift
fill(new_SZ(14).ZONE(:,1),new_SZ(14).ZONE(:,2),[1.0 0.6 0.0],'FaceAlpha',0.3); hold on; %Mweru-1
fill(new_SZ(15).ZONE(:,1),new_SZ(15).ZONE(:,2),[1.0 1.0 1.0],'FaceAlpha',0.7); hold on; %Mweru
fill(new_SZ(16).ZONE(:,1),new_SZ(16).ZONE(:,2),'c','FaceAlpha',0.3); hold on; %Kariba-Okavango-1
fill(new_SZ(17).ZONE(:,1),new_SZ(17).ZONE(:,2),'y','FaceAlpha',0.2); hold on; %Kariba-Okavango-1
fill(new_SZ(19).ZONE(:,1),new_SZ(19).ZONE(:,2),'k','FaceAlpha',0.2); hold on; %Rovuma Basin
fill(new_SZ(20).ZONE(:,1),new_SZ(20).ZONE(:,2),[0.1 0.6 0.0],'FaceAlpha',0.2); hold on; %Eastern and Davie Rift
fill(new_SZ(23).ZONE(:,1),new_SZ(23).ZONE(:,2),[1 1 0],'FaceAlpha',0.2); hold on; %Rovuma Basin
fill(new_SZ(26).ZONE(:,1),new_SZ(26).ZONE(:,2),[0.9 0 0.9],'FaceAlpha',0.3); hold on; %PREPARE SZ_1
fill(new_SZ(27).ZONE(:,1),new_SZ(27).ZONE(:,2),[0.3 0 0.3],'FaceAlpha',0.3); hold on; %PREPARE SZ_2

plot(bg_source{1}(:,2),bg_source{1}(:,3),'.','MarkerEdgeColor',[0.7 0.6 0.2]); hold on
plot(bg_source{2}(:,2),bg_source{2}(:,3),'r.'); hold on
plot(bg_source{3}(:,2),bg_source{3}(:,3),'b.'); hold on
plot(bg_source{4}(:,2),bg_source{4}(:,3),'.','MarkerEdgeColor',[1 0.6 0.0]); hold on
plot(bg_source{5}(:,2),bg_source{5}(:,3),'.','MarkerEdgeColor',[0 0.5 0.8]); hold on
plot(bg_source{6}(:,2),bg_source{6}(:,3),'.','MarkerEdgeColor',[0.7 0.8 0.2]); hold on
plot(bg_source{7}(:,2),bg_source{7}(:,3),'k.'); hold on
plot(bg_source{8}(:,2),bg_source{8}(:,3),'.','MarkerEdgeColor',[0.1 0.6 0.0]); hold on
plot(bg_source{9}(:,2),bg_source{9}(:,3),'.','MarkerEdgeColor',[1 0.4 0]); hold on
plot(bg_source{10}(:,2),bg_source{10}(:,3),'.','MarkerEdgeColor',[0.9 0 0.9]); hold on
plot(bg_source{11}(:,2),bg_source{11}(:,3),'.','MarkerEdgeColor',[0.3 0 0.3]); hold on
plot(buffer_200.X,buffer_200.Y,'r--'); hold on
plot(MapData2(:,2),MapData2(:,1),'k-'); hold on
plot(PSHA_Zone(:,2),PSHA_Zone(:,1),'r','LineWidth',1.5); hold on 
axis equal; axis(Region2); xlabel('Longitude'); ylabel('Latitude')

%% Find MFD of all events in Malawi (i.e. dist<0)
load syncat_bg
rem_event=find(bg_catalog(:,1)>0);
bg_catalog(rem_event,:)=[];

Mo_bg_catalog = 10.^(1.5*bg_catalog(:,4)+9.05);
MoRate_bg_catalog = sum(Mo_bg_catalog)/NumSimu;

