%% Create following SourceCalc Plots:
% Comparison of Zomba Fault and Zomba Graben Fault uncertainity with the MSSM
% Histogram of section, fault, and multifault source magnitudes and recurrence intervals

close all
mydir  = pwd; idcs   = strfind(mydir,'/');
addpath(mydir(1:idcs(end)-1));

%load and sort data

load MSSM_source_calcs
load MSSM_sources

[num_raw,txt_raw,raw_raw]=xlsread('MSSM.xlsx',1);

%clean MSSM source spreadsheet
raw_raw(1:3,:)=[];%Remove header
txt_raw(1:2,:)=[];
raw_raw(height(num_raw)+1:height(raw_raw),:)=[]; %make sure same height as num_raw

%% Obtain slip rate and recurrence interval data for single fault

id = 313; %set for Zomba

%index source and fault name in MSSM spreadsheet
if id >600
        ii=find(id==num_raw(:,83),1);
        f_name=txt_raw(ii+2,6);
    elseif id >300
        ii=find(id==num_raw(:,81),1);
        f_name=txt_raw(ii+2,4);
    else
        ii=find(id==num_raw(:,82),1);
        f_name=txt_raw(ii+2,5);
end

%index source in MSSM_source_calcs
indx = find(id ==MSSM_source_calcs(:,1));

%prob distribution for slip rate
pd_sr1 = makedist('Normal','mu',MSSM_source_calcs(indx,2),'sigma',MSSM_source_calcs(indx,3)); 
pd_r2 = makedist('Lognormal','mu',MSSM_source_calcs(indx,4),'sigma',MSSM_source_calcs(indx,5));
x_sr=linspace(range_sr(indx,1),range_sr(indx,2),1000);
%x_ri=linspace(range_r(indx,1),range_r(indx,2),10e6);
x_ri=linspace(0,range_r(indx,2),5*10e5);

%Select lower, intermediate, and upper estimates from MSSM
if id>600
    x1_sr = [num_raw(ii,38) num_raw(ii,39) num_raw(ii,40)];
    x1_ri = [num_raw(ii,77) num_raw(ii,78) num_raw(ii,79)];
elseif id>300
    x1_sr = [num_raw(ii,41) num_raw(ii,42) num_raw(ii,43)];
    x1_ri = [num_raw(ii,66) num_raw(ii,67) num_raw(ii,68)];
else
    x1_sr = [num_raw(ii,41) num_raw(ii,45) num_raw(ii,46)];
    x1_ri = [num_raw(ii,55) num_raw(ii,56) num_raw(ii,57)];
end

x2_sr = [MSSM_source_calcs(indx,2)-MSSM_source_calcs(indx,3) MSSM_source_calcs(indx,2) MSSM_source_calcs(indx,2)+MSSM_source_calcs(indx,3)];
y1 = [1/3 1/3 1/3];
x2_ri = [exp(MSSM_source_calcs(indx,4)-MSSM_source_calcs(indx,5)) exp(MSSM_source_calcs(indx,4)) exp(MSSM_source_calcs(indx,4)+MSSM_source_calcs(indx,5))];
y2 = [0.16 0.68 0.16];

 
%% Load in data from SMSSD to make comparison with the MSSM

%load in data from SMSSD
[num_SMSSD,txt_SMSSD,raw_SMSSD]=xlsread('SMSSD.xlsx');

zomb_SMSSDindx=find(string(txt_SMSSD(:,3))=='Zomba');
zomb_faults=txt_SMSSD(zomb_SMSSDindx,4); zomb_faults=unique(zomb_faults);
zomb_SMSSDindx=zomb_SMSSDindx-2; %difference between height of txt and num

%find SMSSD slip rates & recurruence intervals
SMSSD_sr=num_SMSSD(zomb_SMSSDindx,31);
SMSSD_srl=num_SMSSD(zomb_SMSSDindx,31)-num_SMSSD(zomb_SMSSDindx,30);%low err
SMSSD_sru=num_SMSSD(zomb_SMSSDindx,32)-num_SMSSD(zomb_SMSSDindx,31);%upp err

SMSSD_ri=num_SMSSD(zomb_SMSSDindx,53);
SMSSD_ril=num_SMSSD(zomb_SMSSDindx,53)-num_SMSSD(zomb_SMSSDindx,52);%low err
SMSSD_riu=num_SMSSD(zomb_SMSSDindx,54)-num_SMSSD(zomb_SMSSDindx,53);%upp err

% Index same faults in MSSM
% Note this will only work for faults in Zomba Graben that are in both MSSM and SMSSD
for i=1:length(zomb_faults)
    if isempty(find(string(zomb_faults(i))==string(txt_raw(:,4)),1))==0
    zomb_MSSMindx(i)=find(string(zomb_faults(i))==string(txt_raw(:,4)),1);
    %for branching faults concatenate -1 so names match
    else
    zomb_MSSMindx(i)=find(string(strcat(zomb_faults(i),'-1'))==string(txt_raw(:,4)),1);  
    end
    
    %match ids with MSSM_source_calcs
    zomba_id(i)=find(num_raw(zomb_MSSMindx(i),1)==MSSM_source_calcs(:,1));
end

%find MSSM slip rates & recurrence intervals
MSSM_sr=MSSM_source_calcs(zomba_id,2); MSSM_sr_err=MSSM_source_calcs(zomba_id,3);

MSSM_ri=exp(MSSM_source_calcs(zomba_id,4));
MSSM_ril=MSSM_ri-exp(MSSM_source_calcs(zomba_id,4)-MSSM_source_calcs(zomba_id,5));
MSSM_riu=exp(MSSM_source_calcs(zomba_id,4)+MSSM_source_calcs(zomba_id,5))-MSSM_ri;



%% Plot combined comparison of single fault and all faults in the Zomba Graben

%Set for comparing particular zomba graben fault
z_id = id; %set for Zomba

zomb_SMSSDindxf=find(string(txt_SMSSD(:,4))==string(txt_raw(find(num_raw(:,1)==z_id,1),4)),1);
x1_sr=[num_SMSSD(zomb_SMSSDindxf,30) num_SMSSD(zomb_SMSSDindxf,31) num_SMSSD(zomb_SMSSDindxf,32)];
x1_ri=[num_SMSSD(zomb_SMSSDindxf,52) num_SMSSD(zomb_SMSSDindxf,53) num_SMSSD(zomb_SMSSDindxf,54)];
%xmax=1+length(zomb_faults)
xmax=1+6;

figure(700);

t=tiledlayout(2,2);
t.TileSpacing = 'compact';
nexttile
%plot slip rate pdf
%subplot(2,2,1)
yyaxis left
ylabel('f(x)'); hold on

y_sr2=pdf(pd_sr1,x_sr);  set(gca,{'ycolor'},{'k'});
plot(x_sr,y_sr2,'k-','linewidth',1.5); hold on

%plot discrete values
yyaxis right
plot(x2_sr,y2,'k*');
plot(x1_sr,y1,'r*');

legend(['Monte Carlo' newline 'simulations' newline '(Normal distribution)'],'Mean slip rate and \pm 1\sigma',...
  ['SMSSD lower,' newline 'intermediate & upper' newline 'slip rate estimates'],'location','northeast','FontSize',10);

ylim ([0 1]); set(gca,{'ycolor'},{'k'}); hold on
xlabel ('slip rate (mm/yr)');  ylabel('Weighting');  hold on
t1=title('(a)','fontsize',16,'fontweight','normal');
set(t1, 'horizontalAlignment', 'left');
set(t1, 'units', 'normalized');
h1 = get(t1, 'position');
set(t1, 'position', [-0.2 h1(2) h1(3)]);

axis square;  set(gca,'fontsize',12); box on

%plot ri pdf
nexttile
%subplot(2,2,2)
yyaxis left
y_ri1=pdf(pd_r2,x_ri);
plot(log10(x_ri),y_ri1,'k-','linewidth',1.5); hold on
ylabel('f(x)'); set(gca,{'ycolor'},{'k'}); 

yyaxis right
plot(log10(x2_ri),y2,'k*');hold on
plot(log10(x1_ri),y1,'r*');hold on

ylim ([0 1]); xlim ([0 8]); set(gca,{'ycolor'},{'k'}); hold on

legend(['Monte Carlo' newline 'simulations (Log' newline 'normal distribution)'],'Mean RI and \pm 1\sigma',...
    ['SMSSD lower,' newline 'intermediate & upper' newline 'RI estimates'],'location','northeast','FontSize',10);
xlabel ('Log10(RI)'); ylabel('Weighting'); hold on


t2=title('(b)','fontsize',16,'fontweight','normal');
set(t2, 'horizontalAlignment', 'left');
set(t2, 'units', 'normalized');
h1 = get(t2, 'position');
set(t2, 'position', [-0.2 h1(2) h1(3)]);


axis square;  set(gca,'fontsize',12); box on


%plot slip rates for all of Zomba Graben
nexttile
%subplot(2,2,3)
errorbar(categorical(zomb_faults),unique(SMSSD_sr),unique(SMSSD_srl),unique(SMSSD_sru),'o',...
    'LineWidth',1.2,'MarkerFaceColor','w','MarkerEdgeColor','r','Color','r');

hold on
errorbar(categorical(zomb_faults),MSSM_sr,MSSM_sr_err,'o',...
    'LineWidth',1.2,'MarkerFaceColor','w','MarkerEdgeColor','k','Color','k');

ylim([0 1]); ylabel('Slip Rate (mm/yr)');
xticklabels(zomb_faults);

t3=title('(c)','fontsize',16,'fontweight','normal');
set(t3, 'horizontalAlignment', 'left');
set(t3, 'units', 'normalized');
h1 = get(t3, 'position');
set(t3, 'position', [-0.2 h1(2) h1(3)]);

set(gca,'fontsize',12,'XTickLabelRotation',(45),'PlotBoxAspectRatio',[1 1 1]);
legend('SMSSD','MSSM','Location','Northwest');

%subplot(2,2,4)
nexttile
%for some reason, flip SMSSD ri to plot in right order (CSF to Zomb)
errorbar(categorical(zomb_faults),flip(unique(SMSSD_ri)),flip(unique(SMSSD_ril)),flip(unique(SMSSD_riu)),'o',...
    'LineWidth',1.2,'MarkerFaceColor','w','MarkerEdgeColor','r','Color','r');

hold on
errorbar(categorical(zomb_faults),MSSM_ri,MSSM_ril,MSSM_riu,'o',...
    'LineWidth',1.2,'MarkerFaceColor','w','MarkerEdgeColor','k','Color','k');

ylim([100 10e6]); ylabel('Recurrence Intervals (Years)')
xticklabels(zomb_faults)

t4=title('(d)','Fontsize',16,'fontweight','normal');
set(t4, 'horizontalAlignment', 'left');
set(t4, 'units', 'normalized');
h1 = get(t4, 'position');
set(t4, 'position', [-0.2 h1(2) h1(3)]);


%Note yaxis set as log scale
set(gca,'yscale','log','fontsize',12,'XTickLabelRotation',(45),'PlotBoxAspectRatio',[1 1 1]);
legend('SMSSD','MSSM','Location','Northwest');


set(gcf,'Position',[681 226 1013 753]);


%% Histogram plots
% Seperately index border and intrabasin faults

% Index unique sources in MSSM spreadsheet by unique magnitudes
% Note this correction does incorrectly remove some sources, so they are
% added back in:

%LIST OF DIFFERENT SECTIONS WITH SAME MAGNITUDES
%Liwawadzi South & Sani Central 
%Malombe South & NB Fault 9a&b North
%Mwanza Thombani, Wamkurumadzi South-2 & Chlingali South
%Sani North, Mbiri-1 South, Kaporo-1 North 
%BMF  Matakataka-1&2, Lweya Central
%Ruo North & Metangula Lake 1&2
%Zomba Chingale Stream & Phirilanyama-2 North
%Thyolo Mbewe & Metangula North 2a
%Central Basin Fault 20 South, Wovwe-2 North, St Mary South
%Bwanje North & South Basin Fault 8 South
%Chirobwe Ncheu-2 Dzonze South, Malombe Chimwalira, South Basin Fault 7b South
%Chingale Step Namitembo, BMF Kasinje 1&2, BMF Lithipe-2
%South Basin Fault 13c South & Central Basin Fault 19 North
%Mwanza Condedezi & South Basin Fault 7 North a&b
%Zomba North & Makanjira Malindi
%Liwawazi North & South Basin Fault 13b Central
%Chirobwe Ncheu Fault Livuzeli 1&2 & Malombe Mpale
%Lipichilli-3 & Mbiri-1 North 
%Thyolo North 1&2 & South Basin Fault 7a South

tmp1=[136,256,60,144,222,232,185,124,142,131,221,250,176,...
    109,172,89,100,205,171,119,165,107,224,170];

%Find sec rows
[fault_sec, ie, ih] = unique(num_raw(:,49));
fault_sec=vertcat(fault_sec,num_raw(tmp1,49));
ie=vertcat(ie,tmp1');

%LIST OF DIFFERENT FAULTS WITH SAME MAGNITUDES
%Bwanwie & Lipichilli North 
%South Fault Basin 6 & Central Basin Fault 7 
%Mbiri-2 & Usisya-North
%South Basin Fault 2a&2b

tmp2=[210,201,190,148];

%Find fault rows 
[fault, ia, ic] = unique(num_raw(:,60));
fault=vertcat(fault,num_raw(tmp2,60));%Correction for equal mag of faults
ia=vertcat(ia,tmp2');

bf_id=zeros(length(fault),1); if_id=zeros(length(fault),1); 

for i=1:length(fault)

if txt_raw(ia(i),28) == "B" && isnan(num_raw(ia(i),42))==0
    bf_id(i)=num_raw(ia(i),1); %index border faults
elseif txt_raw(ia(i),28) == "I" && isnan(num_raw(ia(i),42))==0
    if_id(i)=num_raw(ia(i),1); %index intrabasin faults
    
end
end

bf_id(find(bf_id==0))=[]; if_id(find(if_id==0))=[];
bf_indx = zeros(length(bf_id),1); if_indx = zeros(length(if_id),1);

% Index faults in MSSM_source_calcs
for i=1:length(bf_id)
    bf_indx(i)=find(bf_id(i)==MSSM_source_calcs(:,1));
end

for i=1:length(if_id)
    if_indx(i)=find(if_id(i)==MSSM_source_calcs(:,1));
end

% Index unique sections in MSSM_source_calcs
% Should equal same as num_sec1 in MSSM syncat analysis
sec_indx1=find(MSSM_source_calcs(:,1)<300); tmp_indx=zeros(length(sec_indx1),1);
for i=1:length(sec_indx1)
    tmp_indx(i)=ismember(MSSM_source_calcs(sec_indx1(i),1),num_raw(ie,2));
end

sec_indx1(find(tmp_indx==0))=[]; sec_indx2=zeros(length(sec_indx1),1);

for i=1:length(sec_indx1)
    sec_indx2(i)=find(MSSM_source_calcs(sec_indx1(i),1)==num_sec(:,1));
end

f_indx=vertcat(bf_indx,if_indx);
mf_indx=find(MSSM_source_calcs(:,1)>599);

%% Plot histograms

figure(500);
tiledlayout(1,3,'TileSpacing','Compact');

%Plot slip rates
nexttile
%subplot(1,3,1)
clear h1 h2
h1=histogram(MSSM_source_calcs(bf_indx,2),'BinWidth',0.5,'facecolor',[0.9 0.9 0.9],'facealpha',1); hold on
h2=histogram(MSSM_source_calcs(if_indx,2),'BinWidth',0.05,'facecolor',[0.5 0.5 0.5]); hold on 
legend([h1 h2],{['Border faults, n = ', num2str(length(bf_id))],['Intrarift faults, n = ',num2str(length(if_id))]})
xlim([0 1.6]); ylim([0 60]); ylabel('Frequency'); xlabel('Slip Rate mm/yr'); xticks([0.1 0.2 0.5 1 1.5 2 2.5 3]); xtickangle(-45);...
set(gca,'fontsize',15); axis square; %applyhatch(gcf,'\.'); hold off

%Plot magnitudes 
%Note these are derived from MSSM spreadsheet
nexttile
clear h3 h4 h5
h3=histogram(num_sec(sec_indx2,21),'BinWidth',0.25,'facealpha',0.5); hold on
h4=histogram(unique(num_fault(:,18)),'BinWidth',0.25,'facealpha',0.5); hold on
h5=histogram(num_multi_fault(:,10),'BinWidth',0.25,'facealpha',0.5); hold on
%legend([h3 h4 h5],{['Section ruptures, n = ', num2str(length(sec_indx1))],['Fault ruptures, n = ', num2str(length(f_indx))],['Multi-fault ruptures, n = ',num2str(length(mf_indx))]})
xlabel('Mw'); ylabel('Frequency'); xlim([5.5 8.5]); axis square; set(gca,'fontsize',15); hold off%apply hatch doesn't work %manually do in illustrator%applyhatch(gcf,'./\\.')

%Plot recurrence intervals
nexttile
h6=histogram(log10(exp(MSSM_source_calcs(sec_indx1,4))),'BinWidth',0.5,'facealpha',0.5); hold on
h7=histogram(log10(exp(MSSM_source_calcs(f_indx,4))),'BinWidth',0.5,'facealpha',0.5); hold on
h8=histogram(log10(exp(MSSM_source_calcs(mf_indx,4))),'BinWidth',0.5,'facealpha',0.5); hold on
legend([h6 h7 h8],{['Section ruptures, n = ', num2str(length(sec_indx1))],['Fault ruptures, n = ', num2str(length(f_indx))],['Multi-fault ruptures, n = ',num2str(length(mf_indx))]},'Location','Eastoutside')
xlabel('Recurrence Intervals (years)'); xticks([2 3 4 5]);xticklabels({'10^2','10^3','10^4','10^5'});
ylabel('Frequency'); axis square; set(gca,'fontsize',15);
set(gcf, 'units', 'normalized');
set(gcf, 'Position', [0 0.2044 0.9938 0.6800]);

%% Stats on section, fault, and multifaults

%Slip rates
bfault_sr_stats = [min(MSSM_source_calcs(bf_indx,2)) mean(MSSM_source_calcs(bf_indx,2)) max(MSSM_source_calcs(bf_indx,2))];
ifault_sr_stats = [min(MSSM_source_calcs(if_indx,2)) mean(MSSM_source_calcs(if_indx,2)) max(MSSM_source_calcs(if_indx,2))];

%Magnitudes
sec_mag_stats = [min(num_sec(sec_indx2,21)) mean(num_sec(sec_indx2,21)) max(num_sec(sec_indx2,21))];
fault_mag_stats = [min(unique(num_fault(:,18))) mean(unique(num_fault(:,18))) max(unique(num_fault(:,18)))];
mfault_mag_stats =  [min(num_multi_fault(:,10)) mean(num_multi_fault(:,10)) max(num_multi_fault(:,10))];

%Recurrence Intervals
sec_ri_stats = [min(exp(MSSM_source_calcs(sec_indx1,4))) mean(exp(MSSM_source_calcs(sec_indx1,4))) max(exp(MSSM_source_calcs(sec_indx1,4)))];
fault_ri_stats = [min(exp(MSSM_source_calcs(f_indx,4))) mean(exp(MSSM_source_calcs(f_indx,4))) max(exp(MSSM_source_calcs(f_indx,4)))];
mfault_ri_stats =  [min(exp(MSSM_source_calcs(mf_indx,4))) mean(exp(MSSM_source_calcs(mf_indx,4))) max(exp(MSSM_source_calcs(mf_indx,4)))];

%Ruptures >Mw 7.5
tmp_mag_store=vertcat(num_sec(:,21),num_fault(:,18),num_multi_fault(:,10));
large_e= length(find(tmp_mag_store>7.5));



%% ARCHIVE

%Plot slip rate and R comparison for just single fault
figure(101);

%plot slip rate pdf
subplot(1,2,1)
yyaxis left
ylabel('f(x)'); hold on

y_sr2=pdf(pd_sr1,x_sr); 
plot(x_sr,y_sr2,'k-','linewidth',1.5); hold on

%plot discrete values
yyaxis right
plot(x2_sr,y2,'k*');
plot(x1_sr,y1,'r*');

legend(['Monte Carlo simulation results' newline '(Normal distribution)'],'Mean slip rate and \pm 1\sigma',...
  ['MSSM lower, intermediate,' newline ' & upper slip rate estimates'],'location','northeast');

ylim ([0 1]); hold on
xlabel ('slip rate (mm/yr)');  ylabel('Weighting'); title([string(f_name) 'Slip Rate']); subtitle('10000 simulations'); hold on
axis square;  set(gca,'fontsize',13); box on

%plot ri pdf
subplot(1,2,2)
yyaxis left
y_ri1=pdf(pd_r2,x_ri);
plot(log10(x_ri),y_ri1,'k-','linewidth',1.5); hold on
ylabel('f(x)'); 

yyaxis right
plot(log10(x2_ri),y2,'k*');hold on
plot(log10(x1_ri),y1,'r*');hold on

ylim ([0 1]); hold on

legend(['Monte Carlo simulation results' newline '(Log normal distribution)'],'Mean RI and \pm 1\sigma',...
    ['MSSM lower, intermediate,' newline ' & upper RI estimates'],'location','northeast');
xlabel ('Log10(RI)'); ylabel('Weighting'); title([string(f_name) 'Recurrence Interval']); subtitle('10000 simulations');  hold on

axis square;  set(gca,'fontsize',13); box on

set(gcf,'position',[788 275 1133 519]) 


%% Plot slip rate comparisons within the Zomba Graben

figure(600);

%plot slip rates
subplot(1,2,1)
errorbar(categorical(zomb_faults),unique(SMSSD_sr),unique(SMSSD_srl),unique(SMSSD_sru),'o',...
    'LineWidth',1.2,'MarkerFaceColor','w','MarkerEdgeColor','k','Color','k');

hold on
errorbar(categorical(zomb_faults),MSSM_sr,MSSM_sr_err,'o',...
    'LineWidth',1.2,'MarkerFaceColor','w','MarkerEdgeColor','r','Color','r');

%xlim([0 1+length(zomb_faults)]); xticks([1:1:length(zomb_faults)]); xticklabels(zomb_faults);
ylim([0 1]); ylabel('Slip Rate (mm/yr)')
t1 = title('(a)');
set(t1, 'horizontalAlignment', 'left');
set(t1, 'units', 'normalized');
h1 = get(t1, 'position');
set(t1, 'position', [-0.2 h1(2) h1(3)]);
set(gca,'fontsize',12,'XTickLabelRotation',(45),'PlotBoxAspectRatio',[1 1 1]);
legend('SMSSD','MSSM','Location','Northwest');

%plot recurrence interval
subplot(1,2,2)
%for some reason, flip SMSSD ri to plot in right order (CSF to Zomb)
errorbar(categorical(zomb_faults),flip(unique(SMSSD_ri)),flip(unique(SMSSD_ril)),flip(unique(SMSSD_riu)),'o',...
    'LineWidth',1.2,'MarkerFaceColor','w','MarkerEdgeColor','k','Color','k');

hold on
errorbar(categorical(zomb_faults),MSSM_ri,MSSM_ril,MSSM_riu,'o',...
    'LineWidth',1.2,'MarkerFaceColor','w','MarkerEdgeColor','r','Color','r');

%xlim([0 1+length(zomb_faults)]); 
%xticks([1:1:length(zomb_faults)]); xticklabels(zomb_faults), 
ylabel('Recurrence Intervals (Years)'); ylim([100 10e6]);
t2 = title('(b)');
set(t2, 'horizontalAlignment', 'left');
set(t2, 'units', 'normalized');
h1 = get(t2, 'position');
set(t2, 'position', [-0.2 h1(2) h1(3)]);


%Note yaxis set as log scale
set(gca,'yscale','log','fontsize',12,'XTickLabelRotation',(45),'PlotBoxAspectRatio',[1 1 1]);
legend('SMSSD','MSSM','Location','Northwest');

set(gcf,'Position',[681 566 792 413]);