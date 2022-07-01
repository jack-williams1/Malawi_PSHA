%% Plot flexural profiles through basins including topography and faults

%PLOTTING TOPOGRAPHIC PROFILES REQUIRES:
%1.) TOPOTOOLBOX (https://au.mathworks.com/matlabcentral/fileexchange/50124-topotoolbox)
%2.) TanDEM-X DATA FOR MALAWI, AVAILABLE BY REQUEST FROM THE GERMAN
%AEROSPACE CENTRE (DLR) AT https://tandemx-science.dlr.de/cgi-bin/wcm.pl?page=TDM-Proposal-Submission-Procedure
%3.) SRTM 30 M DEM FOR MALAWI 
%4.) MATLAB FUNCTION BOUNDED LINE (https://au.mathworks.com/help/matlab/ref/exist.html)

%load data from flexure calculations
load hangingwallflexcalcs;

num_data=table2array(num_table(:,2:width(num_table)));
Wo = [num_data(:,12)-num_data(:,13) num_data(:,12) num_data(:,12)+num_data(:,13)];

% Read in data with projected distance between border and intrarift fauls
dist_2bf=readtable('HangingWallFlexureInputs.xlsx','Sheet','dist2bf');


%% Add topographic swaths across basins

%EDIT PATH AS NECESSARY TO TOPOTOOLBOX AND TANDEMX FOLDERS
addpath(genpath('/Users/wilja48p/Documents/MATLAB')); %for desktop
%addpath(genpath('/Users/jackwilliams/Documents/topotoolbox')); %for macbook

%addpath('/Users/jackwilliams/Documents/Cardiff/TerrainGISFiles_Malawi/TanDEM-X');
%addpath('/Users/jackwilliams/Documents/Cardiff/TerrainGISFiles_Malawi/GMTSAR_SRTM');
addpath('CrossSectionProfiles')

%Read in cross section profiles in qgis
nbasin_section=shaperead('NorthBasin.shp'); cbasin_section=shaperead('CentralBasin.shp'); 
sbasin_section=shaperead('SouthBasin.shp'); makanjira_section=shaperead('Makanjira.shp');
zomba_section=shaperead('Zomba.shp'); lowershire_section=shaperead('LowerShire.shp');

s_malawi_basins=[makanjira_section zomba_section lowershire_section];
l_malawi_basins=[nbasin_section cbasin_section sbasin_section]; %add more basins as appropriate

%{
% Create swaths, takes long time to read in DEM's so only run once
tanDEMn = GRIDobj('merged_north.tif'); %Read in TanDEM-X DEM


tanDEMn_utm=reproject2utm(tanDEMn,12); %Project DEM into utm

tanDEMs = GRIDobj('merged_south.tif'); %Read in TanDEM-X DEM
tanDEMs_utm=reproject2utm(tanDEMs,12); %Project DEM into utm

% Use SRTM for North Basin as TanDEM-X coverage insufficient
srtm_dem=GRIDobj('Malawi_GMTSAR_Merged.tif');
srtm_utm=reproject2utm(srtm_dem,30); %Project DEM into utm

for i=1:length(l_malawi_basins)
    clrX=isnan(l_malawi_basins(i).X);  clrY=isnan(l_malawi_basins(i).Y);
    l_malawi_basins(i).X(clrX)=[]; l_malawi_basins(i).Y(clrY)=[];
    
    clear clrX clrY
end

%flip central basin so all profiles are from west to east
l_malawi_basins(2).X=flip(l_malawi_basins(2).X); l_malawi_basins(2).Y=flip(l_malawi_basins(2).Y);

%Create swathe along basin section with 10 km 
%Note can't create lake basin swaths in for loop as they require different DEMS

nbasinSW=SWATHobj(srtm_utm,l_malawi_basins(1).X,l_malawi_basins(1).Y,'width',10000);
cbasinSW=SWATHobj(tanDEMn_utm,l_malawi_basins(2).X,l_malawi_basins(2).Y,'width',10000);
sbasinSW=SWATHobj(tanDEMs_utm,l_malawi_basins(3).X,l_malawi_basins(3).Y,'width',10000);

lbasinSW=[nbasinSW cbasinSW sbasinSW];
save('l_basinSW.mat','lbasinSW');

%% For south Malawi Basins

for i=1:length(s_malawi_basins)
    clrX=isnan(s_malawi_basins(i).X);  clrY=isnan(s_malawi_basins(i).Y);
    s_malawi_basins(i).X(clrX)=[]; s_malawi_basins(i).Y(clrY)=[];
    
    clear clrX clrY
    s_malawi_swath(i)=SWATHobj(tanDEMs_utm, s_malawi_basins(i).X,s_malawi_basins(i).Y,'width',10000);
    save('s_malawi_swaths','s_malawi_swath')
end
%}

%load topo swaths
load l_basinSW.mat; load s_malawi_swaths

%% Clean lake profiles of noisy elevations around lake and shift profiles

l_elev=460;

lbf_shift=[15000, 10000, 24500];%distance in m between end of section and border fault trace
bf_name={'Livingstone','Usisya','SB 5-13'};
dip=[60 53 53]; lbsin_indx=[3,4,5]; f_idx={}; fy_dist={};


for i=1:length(lbasinSW)

lake_p=find(lbasinSW(i).Z<l_elev); lbasinSW(i).Z(lake_p)=l_elev;
f_idx{i}=find(string(table2cell(num_table(lbsin_indx(i),1)))==string(table2cell(dist_2bf(:,2))));

if i~=2 %shift is +ve for Central Basin (profile starts at 0)
 
    dist_shift(i)=lbf_shift(i)+Wo(i+2,2)*cos(dip(i)*pi/180);
    prof_shift{i}=distance+dist_shift(i);
   
    %position of intrarift faults
    fx_dist{i}=dist_shift(i)+table2array(dist_2bf(f_idx{i},4))*1000;

else  %shift is -ve for central basins (profile starts at max distx)
      
    dist_shift(i)=max(lbasinSW(i).distx)-lbf_shift(i)-Wo(i+2,2)*cos(dip(i)*pi/180);
    prof_shift{i}=dist_shift(i)-distance;

    %position of intrarift faults
    fx_dist{i}=dist_shift(i)-table2array(dist_2bf(f_idx{i},4))*1000;
    %position of border fault
    x_tmp=max(lbasinSW(i).distx)-lbf_shift(i);
end


%index poisiton of fault along topo profile
%note need to flip topo profile
for j=1:height(fx_dist{i})
    [~,dist_indx(j)]=min(abs(lbasinSW(i).distx-fx_dist{i}(j)));
end

z_mean = nanmean(lbasinSW(i).Z,1)';
fy_dist{i}=z_mean(dist_indx);

clear lake_p dist_indx z_mean
end


%% Clean profiles and sort data for southern basins

sbf_shift=[16100,22800,21000];%distance in m between end of section and border fault trace
sbf_name={'Chriobwe-Ncheu','Zomba','Thyolo'};
sbsin_indx=[6,8,9];

for i=1:length(s_malawi_swath)
    
    f_idx{i}=find(string(table2cell(num_table(sbsin_indx(i),1)))==string(table2cell(dist_2bf(:,2))));
   
    
    if i==1  %only correct elevations for Makanjira across Lake Malombe
    lake_p=find(s_malawi_swath(i).Z<460); s_malawi_swath(i).Z(lake_p)=460;    
     
    %shift is +ve for Makanjira (profile starts at 0)
    dist_shift_s(i)=sbf_shift(i)+Wo(sbsin_indx(i),2)*cos(53*pi/180);
    prof_shift_s{i}=distance+dist_shift_s(i);
    
    %sort out shift for plotting border faults
    [~,bfdist_indx(i)]=min(abs(s_malawi_swath(i).distx-sbf_shift(i)));
    x_tmp_s=dist_shift_s(i)+90000;
    [~,bfdist_indx2]=min(abs(s_malawi_swath(i).distx-x_tmp)); 
   
    %position of intrarift faults
    fx_dist_s{i}=table2array(dist_2bf(f_idx{i},4))*1000+dist_shift_s(i);
   
    
    else
    %shift is -ve for other basins (profile starts at max distx)
    [~,bfdist_indx(i)]=min(abs(s_malawi_swath(i).distx-(max(s_malawi_swath(i).distx)-sbf_shift(i))));
    dist_shift_s(i)=max(s_malawi_swath(i).distx)-sbf_shift(i)-Wo(sbsin_indx(i),2)*cos(53*pi/180);
    prof_shift_s{i}=dist_shift_s(i)-distance;
    %position of intrarift faults
    fx_dist_s{i}=dist_shift_s(i)-table2array(dist_2bf(f_idx{i},4))*1000;
     
    end
    
    %index poisiton of fault along topo profile
    for j=1:height(fx_dist_s{i})
        [~,dist_indx_s(j)]=min(abs(s_malawi_swath(i).distx-fx_dist_s{i}(j)));
    end
    
    z_mean = nanmean(s_malawi_swath(i).Z,1)';
    sbfy_dist(i)=z_mean(bfdist_indx(i));
    fy_dist_s{i}=z_mean(dist_indx_s);
    
    if i==1 %indx position of makanjira
       y_tmp_s=z_mean(bfdist_indx2);
    end
    
    clear z_mean dist_indx_s
    
%correct  xaxis s-fx_dist{i}(j)));
end


%% Plot flexural profile through lake basins with topography and faults

%plot with 3x vertical exagg
vegg=3;

rotation=[-99,-81,-99]; ax3_pos=zeros(length(lbasinSW),4);
tmpshift=[6000,4000,5500]; yplotmin=-6000;

clear ax3 ax4

figure(1001);

for i=1:length(l_malawi_basins)

subplot(length(l_malawi_basins),1,i);
    
%Plot topographic cross section
plotdz(lbasinSW(i),'boundedline',true);hold on

%Plot projection of border fault
if i~=2 
line([lbf_shift(i)-yplotmin*cos(dip(i)*pi/180),lbf_shift(i)],[yplotmin,l_elev],'linestyle','--','Color','k');hold on
bf_label=text(lbf_shift(i)-2000,yplotmin+tmpshift(i),bf_name{i},'FontSize',9);
else
line([x_tmp+yplotmin*cos(dip(i)*pi/180),x_tmp],[yplotmin,l_elev],'linestyle','--','Color','k');hold on 
bf_label=text(dist_shift(i)+3500,yplotmin+tmpshift(i),bf_name{i},'FontSize',9);
end

set(bf_label,'Rotation',rotation(i));

if i==3
line([11500-yplotmin*cos(dip(i)*pi/180),11500],[yplotmin,l_elev],'linestyle','--','Color','k');hold on 
bf_label=text(11500-2000,yplotmin+tmpshift(i),'Metangula','FontSize',9); 
set(bf_label,'Rotation',rotation(3));
end


%fault projection
plot1=plot(fx_dist{i},fy_dist{i},'k*');hold on

line1=line(prof_shift{i},squeeze(profile(i+2,2,:))+l_elev,'linewidth',1,'color','b'); %Plot median hanging wall profile
line2=line(prof_shift{i},squeeze(profile(i+2,1,:))+l_elev,'LineStyle','--','linewidth',1,'color','b'); %Plot min hanging wall profile
line3=line(prof_shift{i},squeeze(profile(i+2,3,:))+l_elev,'LineStyle','--','linewidth',1,'color','b'); %Plot maxn hanging wall profile


if i~=2 
xticks([dist_shift(i)-20000:20000:dist_shift(i)+60000]); xticklabels({'-20','0','20','40','60'});
else
xticks([dist_shift(i)-60000:20000:dist_shift(i)]); xticklabels({'60','40','20','0'});
end

yticks([yplotmin:2000:2000]); yticklabels({'-6','-4','-2','0','2'});

ax3(i)=gca;  set(ax3(i),'FontSize',9); set(ax3(i),'XDir','reverse','Units','pixels');  
xlim(ax3(i),[0 max(lbasinSW(i).distx)]); ylim(ax3(i),[yplotmin 3000]); 
pbaspect(ax3(i),[max(lbasinSW(i).distx)/(9000*vegg),1,1]);
%daspect(ax3(i),[max(lbasinSW(i).distx)/4*9000 9000  1]);
ylabel(ax3(i),('Elevation (km)')); xlabel(ax3(i),('Distance from border fault (km)'))


ax3(i).Position(4)=79.2610; %weird fix to make sure all figures all same height


ax3_pos(i,:)=ax3(i).Position;

ax4(i) = axes('Units','pixels','Position',ax3_pos(i,:),'XAxisLocation','top',...
           'YaxisLocation','right','Color','none','xtick',[]);%Make new set of axes for strain plot

%daspect(ax4(i),[max(lbasinSW(i).distx)/4*10000 1  1]); 
 
xlim(ax4(i),[0 max(lbasinSW(i).distx)]);  ylim(ax4(i),[0 10]);   

%Plot strain profiles  
line4=line(prof_shift{i},squeeze(strain_pc(i+2,2,:)),'color','red'); %Plot median hanging wall profile       
line5=line(prof_shift{i},squeeze(strain_pc(i+2,1,:)),'LineStyle','--','color','red'); hold on %Plot min hanging wall profile
line6=line(prof_shift{i},squeeze(strain_pc(i+2,3,:)),'LineStyle','--','color','red'); hold on %Plot max hanging wall profile   
ylabel (ax4(i), 'Strain (%)'); set(ax4(i),'FontSize',9); set(ax4(i),'XDir','reverse');   
pbaspect(ax4(i),[max(lbasinSW(i).distx)/(9000*vegg),1,1]);
   
end

l=legend([plot1 line1 line4],{['Intrarift fault' char(10) 'projection'],['Hanging wall' char(10) 'profile'],'Strain'},...
    'Units','pixels','Location','Eastoutside','box','off','FontSize',9);
l.Position=[515 62.5 93 60];%fix position of legend outside box so doesn't resize plots

set(gcf,'Position',[681 559 621 420],'Units','pixels');
%set(gcf,'Position',[681 425 916 554]);

%% Plot flexural profile through south basins with topography and faults

ax3_pos=zeros(length(s_malawi_swath),4);
ax3_outerpos=zeros(length(s_malawi_swath),4);

figure(1003);

clear ax3 ax4

%plot with 6x vertical exagg
yplotmin_s=[-1000, -1000, -2000]; vegg=6;

for i=1:length(s_malawi_swath)
    
    
    subplot(length(s_malawi_swath),1,i)
    
    %Plot topographic cross section
    plotdz(s_malawi_swath(i),'boundedline',true);hold on
    
    plot1=plot(fx_dist_s{i},fy_dist_s{i},'k*');hold on
    
   
    if i==1 %correction for makanjira
    
    %plot border fault traces (need to do for CNF and Makanjira)    
    line([sbf_shift(i)-yplotmin_s*cos(53*pi/180),sbf_shift(i)],[yplotmin_s,sbfy_dist(i)],'linestyle','--','Color','k');hold on 
    bf_label1=text(dist_shift(i)-4000,-600,['Chirobwe' char(10) 'Ncheu'],'FontSize',8);
    set(bf_label1,'Rotation',92);
    
    
    line([x_tmp_s+yplotmin_s*cos(53*pi/180),x_tmp_s],[yplotmin_s,y_tmp_s],'linestyle','--','Color','k');hold on 
    bf_label2=text(x_tmp_s+2000,500,'Makanjira','FontSize',8);
    set(bf_label2,'Rotation',268);
    
   	line1=line(prof_shift_s{i}(1:90),mak_totProfile(2,:)+l_elev,'linewidth',1,'color','b'); %Plot median hanging wall profile
    line2=line(prof_shift_s{i}(1:90),mak_totProfile(1,:)+l_elev,'LineStyle','--','linewidth',1,'color','b'); %Plot min hanging wall profile
    line3=line(prof_shift_s{i}(1:90),mak_totProfile(3,:)+l_elev,'LineStyle','--','linewidth',1,'color','b'); %Plot maxn hanging wall profile
    
    xticks([dist_shift_s(i):20000:dist_shift_s(i)+100000]); xticklabels({'0','20','40','60','80','100'});
    
    elseif i==2 %Plot Zomba
    
    line([max(s_malawi_swath(i).distx)-sbf_shift(i)+yplotmin_s*cos(53*pi/180),max(s_malawi_swath(i).distx)-sbf_shift(i)],[yplotmin_s,sbfy_dist(i)],'linestyle','--','Color','k');hold on     
    bf_label3=text(dist_shift_s(i)+2000,400,sbf_name(i),'FontSize',8);
    set(bf_label3,'Rotation',268);  
    
    
    line1=line(prof_shift_s{i},squeeze(profile(sbsin_indx(i),2,:))+l_elev,'linewidth',1,'color','b'); %Plot median hanging wall profile
    line2=line(prof_shift_s{i},squeeze(profile(sbsin_indx(i),1,:))+l_elev,'LineStyle','--','linewidth',1,'color','b'); %Plot min hanging wall profile
    line3=line(prof_shift_s{i},squeeze(profile(sbsin_indx(i),3,:))+l_elev,'LineStyle','--','linewidth',1,'color','b'); %Plot maxn hanging wall profile
  
    xticks([dist_shift_s(i)-80000:20000:dist_shift_s(i)]); xticklabels({'80','60','40','20','0'});
    
    else %Lower Shire at ~100 m elevation, needs differnet elevation correction
    
    line([max(s_malawi_swath(i).distx)-sbf_shift(i)+yplotmin_s(i)*cos(53*pi/180),max(s_malawi_swath(i).distx)-sbf_shift(i)],[yplotmin,sbfy_dist(i)],'linestyle','--','Color','k');hold on  
    bf_label4=text(dist_shift_s(i)+2000,0,sbf_name(i),'FontSize',8);
    set(bf_label4,'Rotation',268);  
        
    line1=line(prof_shift_s{i},squeeze(profile(sbsin_indx(i),2,:))+100,'linewidth',1,'color','b'); %Plot median hanging wall profile
    line2=line(prof_shift_s{i},squeeze(profile(sbsin_indx(i),1,:))+100,'LineStyle','--','linewidth',1,'color','b'); %Plot min hanging wall profile
    line3=line(prof_shift_s{i},squeeze(profile(sbsin_indx(i),3,:))+100,'LineStyle','--','linewidth',1,'color','b'); %Plot maxn hanging wall profile    
        
    xticks([dist_shift_s(i)-60000:20000:dist_shift_s(i)]); xticklabels({'60','40','20','0'});
        
    end
    ax3(i)=gca;
    
    xlim(ax3(i),[0 max(s_malawi_swath(i).distx)]); ylim(ax3(i),[yplotmin_s(i) 2000]);
    ylabel(ax3(i),('Elevation (km)')); xlabel(ax3(i),('Distance from border fault (km)'));
    if i~=3
        yticks([-1000:1000:2000]); yticklabels({'-1','0','1','2'});
    else
        yticks([-2000:1000:2000]); yticklabels({'-2','-1','0','1','2'});
    end
    set(ax3(i),'FontSize',9,'Units','pixels'); yrange=(2000-yplotmin_s(i));
    %daspect(ax3(i),[1 6/(max(s_malawi_swath(1).distx)/yrange) 1]); 
    pbaspect(ax3(i),[max(s_malawi_swath(i).distx)/(yrange*vegg),1,1]);
    %make new axes for strain plots
    
    ax3_pos(i,:)=ax3(i).Position; 
    ax4(i) = axes('Units','pixels','Position',ax3_pos(i,:),'XAxisLocation','top',...
           'YaxisLocation','right','Color','none','xtick',[]);%'PositionConstraint','innerposition');%Make new set of axes for strain plot
    xlim(ax4(i),[0 max(s_malawi_swath(i).distx)]); ylim(ax4(i),[0 4]);
    %daspect(ax4(i),[1 6/(max(s_malawi_swath(1).distx)/4) 1]);
    pbaspect(ax4(i),[max(s_malawi_swath(i).distx)/(yrange*vegg),1,1]);
     if i==1  
                  
    line4=line(ax4(i),prof_shift_s{i}(1:90),mak_totStrain_pc(2,:),'color','red');
    line5=line(ax4(i),prof_shift_s{i}(1:90),mak_totStrain_pc(1,:),'LineStyle','--','color','red'); %Plot min hanging wall profile
    line6=line(ax4(i),prof_shift_s{i}(1:90),mak_totStrain_pc(3,:),'LineStyle','--','color','red'); %Plot maxn hanging wall profile 
     else
         
    line4=line(ax4(i),prof_shift_s{i},squeeze(strain_pc(sbsin_indx(i),2,:)),'color','red'); %Plot median hanging wall profile      
    line5=line(ax4(i),prof_shift_s{i},squeeze(strain_pc(sbsin_indx(i),1,:)),'LineStyle','--','color','red'); hold on %Plot min hanging wall profile
    line6=line(ax4(i),prof_shift_s{i},squeeze(strain_pc(sbsin_indx(i),3,:)),'LineStyle','--','color','red'); hold on %Plot max hanging wall profile   

     end
     
    set(ax4(i),'FontSize',9); ylabel(ax4(i),('Strain (%)')); 

end

l=legend([plot1(1) line1 line4],{['Intrarift fault' char(10) 'projection'],['Hanging wall' char(10) 'profile'],'Strain'},...
    'Units','pixels','Location','Eastoutside','box','off','FontSize',9);
l.Position=[515 42.5 93 60];%fix position of legend outside box so doesn't resize plots

set(gcf,'Position',[681 559 621 420]);