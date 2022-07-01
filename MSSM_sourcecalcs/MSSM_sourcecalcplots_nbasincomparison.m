%% Create various types of plots to depict epistemic uncertainty

close all

%load and sort data

load MSSM_source_calcs

addpath('nbasin_faultdata');

nbasin_faults_sr=readtable('nb_slipratetable_id.csv');
%load spreadsheet that cross references MSSM and Shillington et al 2020 fault names
nbasin_faults_xref=readtable('nbasin_faults_idxref.xlsx');

nbasin_faults=table2array(nbasin_faults_sr);
%load slip rates from Shillington et al 2020 calculated through Frankel and Zechar (2009) matlab scripts
load('v_id.mat');
mydir  = pwd; idcs   = strfind(mydir,'/');
addpath(mydir(1:idcs(end)-1));
addpath([mydir(1:idcs(end)-1) '/misc_functions']);
%% Plot slip rate and recurrence interval variation as a normal distribution
%plots with intital MSSM estimates and mean and +/-1 sigma

id = 414; %id of MSSM source to plot slip rate and recurrence interval of
           %Make sure id is one of the nb faults
    
    
%index source in MSSM_source_calcs
indx = find(id ==MSSM_source_calcs(:,1));

%index source in nbasin_faults
l_sr_indx=find(id==cell2mat(slip_rate_nb.id));
%index slip rate from Shillington 2020
l_sr_indx2=nbasin_faults_xref.Fault_Shillington2020_(find(id==nbasin_faults_xref.MSSM_id));

%make normal distribution from simulations
pd_sr1 = makedist('Normal','mu',mean(slip_rate_nb.sr{l_sr_indx}),'sigma',std(slip_rate_nb.sr{l_sr_indx}));
x_sr1=linspace(min(slip_rate_nb.sr{l_sr_indx}),max(slip_rate_nb.sr{l_sr_indx}),1000);

%prob distribution for slip rate form offset 75 Ka (as used in MSSM_source_calcs)
pd_sr2 = makedist('Normal','mu',MSSM_source_calcs(indx,2),'sigma',MSSM_source_calcs(indx,3)); 
x_sr2=linspace(min(v{l_sr_indx2}),max(v{l_sr_indx2}),1000); %use x_sr values from Zechar and Frankel (2009)

%Obtain 1 sigma errorbars for all faults
% y axis: Shillington Values
yerr = nbasin_faults(:,4);

% x axis: MSSM values-Simulation Error
%load in values from slip_rate_nb.sr generated during simulations
for nn=1:height(nbasin_faults)
    
    nbasin_faults_index(nn)=find(nbasin_faults(nn,1)==cell2mat(slip_rate_nb.id));
    simu_mean(nn)=mean(slip_rate_nb.sr{nbasin_faults_index(nn)});%mean value from simulations
    xerr(nn)=std(slip_rate_nb.sr{nbasin_faults_index(nn)});%1 standard dev from simulations
  
end

%% Plot figure with comparison of (a) slip rate pdf for single fault
% and (b) comparison of all faults in scatter plot
figure(101);

subplot(1,2,1)
ylabel('f(x)'); hold on

%plot MSSM pdf
y_sr1=pdf(pd_sr1,x_sr1); 
plot(x_sr1,y_sr1,'k-','linewidth',1.5); hold on

%plot reflector slip rate pdf
y_sr2=pdf(pd_sr2,x_sr2);
plot(x_sr2,y_sr2,'linewidth',1.5); hold on

legend(['Slip rate: Monte Carlo simulation' newline '(Normal distribution)'],['Slip rate: offset 75Ka reflector' newline '(Normal distribution)'],...
  'location','northeast');

%ylim ([0 1]); hold on
xlabel ('slip rate (mm/yr)'); %title([string(f_name) 'Slip Rate']); subtitle('10000 simulations'); hold on
axis square;  set(gca,'fontsize',13); box on

subplot(1,2,2)

scatter(simu_mean,nbasin_faults(:,3)); hold on;
errorbar(simu_mean,nbasin_faults(:,3),yerr,yerr,xerr,xerr,'o','LineWidth',1,'MarkerFaceColor','w','MarkerEdgeColor','k','Color','k'); hold on
%errorbar(nbasin_faults(:,2),MSSM_source_calcs(nbasin_faults_index,2),yneg,ypos,xerr,xerr,'o','LineWidth',1,'MarkerFaceColor','w','MarkerEdgeColor','r','Color','r'); hold on;
plot([0 1],[0 1],'k--'); annotation('textarrow',[0.8 0.73],[0.35 0.5],'String','y = x ','FontSize',13); set(gca,'FontSize',13); hold on; %Plot line of y=x 
xlim([0 1]); ylim([0 1]); xlabel('MSSM Slip Rate (mm/yr)'); ylabel('Shillington et al (2020) Slip Rate (mm/yr)'); axis square; box on

set(gcf,'position',[788 275 1133 519])  

%% Plot figure comparing slip rate estimate pdfs and calculate overlap coefficient
figure(1001);
indx=0;

t=tiledlayout(3,4);
t.TileSpacing = 'compact';
nsim=600; %number of monte carlo simulations

for ii=1:height(nbasin_faults)
    
    %only index first case to avoid MSSM duplications (e.g. nbasin fault 4)
    if ii==find(nbasin_faults(ii,3)==nbasin_faults(:,3),1)
    
    indx=indx+1;  
    
    %find title of fault
    fault_name=nbasin_faults_xref.Fault_Scholz2020_MSSM_(find(nbasin_faults(ii,1)==nbasin_faults_xref.MSSM_id));
    
    %10000 slip rates from thset(gcf,'position',[788 275 1133 519])  
    simu_samp=normrnd(simu_mean(ii),xerr(ii),nsim,1);
    sr_samp=normrnd(nbasin_faults(ii,3),yerr(ii),nsim,1);
    
    %enforce positivity constraint
    simu_samp(find(simu_samp<0))=0.001;  sr_samp(find(sr_samp<0))=0.001; 
    
    [h_tmp(indx),pval(indx),test(indx),df(indx)]= chi2test2(simu_samp,sr_samp);
    
    %Option to run  two sample t-test with unequal variance
    %[h(indx) pval(indx)] = ttest2(sr_samp,simu_samp,'Vartype','unequal');
    %Option to run ztest, but note this assumes each distribution has the same variance
    %[h(indx) pval(ii) ci(:,ii) zval(ii)] = ztest(simu_samp,nbasin_faults(ii,3),yerr(ii));
    
    nexttile
    [ovlr(indx),h1(indx), h2(indx), h3(indx)]=calc_overlap_twonormal(xerr(ii),yerr(ii),simu_mean(ii),nbasin_faults(ii,3),0,1.2,0.01,...
    pval(indx),fault_name,yerr(ii));
    
    end
   
    hold off
end

%add legend
Lgnd=legend([h1(indx) h2(indx) h3(indx)],{'Systems-based','Offset Reflector','Overlap'},'FontSize',13);
Lgnd.Position(1) = 0.75;
Lgnd.Position(2) = 0.24;


set(gcf,'Position', [681 84 1083 895]);
%%

indx=0; nsim2=500;

t=tiledlayout(3,4);
t.TileSpacing = 'compact';
nsim=[100,1000,10000]; %number of monte carlo simulationd
h=zeros(nsim2,11);

for kk=1:length(nsim)
 
for jj=1:nsim2   
    
for ii=1:height(nbasin_faults)
    
    %only index first case to avoid MSSM duplications (e.g. nbasin fault 4)
    if ii==find(nbasin_faults(ii,3)==nbasin_faults(:,3),1)
    
    indx=indx+1;  
    
    %find title of fault
    fault_name=nbasin_faults_xref.Fault_Scholz2020_MSSM_(find(nbasin_faults(ii,1)==nbasin_faults_xref.MSSM_id));
    
    %10000 slip rates from thset(gcf,'position',[788 275 1133 519])  
    simu_samp=normrnd(simu_mean(ii),xerr(ii),nsim(kk),1);
    sr_samp=normrnd(nbasin_faults(ii,3),yerr(ii),nsim(kk),1);

    [h(jj,indx) pval] = ttest2(sr_samp,simu_samp,'Vartype','unequal');
  
    end
    hold off
end %end ii loop for each fault

  indx=0;
  
end %end jj loop for nsim2

%For each kk loop, find number of tests that 
%cannot reject the null hypothesis
for gg=1:11
    
   h_count(gg,kk)=length(find(h(:,gg)<1));
   
end

end %end kk loop for nsim

%% Plot results


indx=0;

%Get fault names in single array
for ii=1:height(nbasin_faults)
    
    %only index first case to avoid MSSM duplications (e.g. nbasin fault 4)
    if ii==find(nbasin_faults(ii,3)==nbasin_faults(:,3),1)
    
    indx=indx+1;  
    
    %find title of fault
    fault_name(indx)=nbasin_faults_xref.Fault_Scholz2020_MSSM_(find(nbasin_faults(ii,1)==nbasin_faults_xref.MSSM_id));
    end
end

figure(2);
bar(h_count); hold on;
legend({'100 simulations','1000 simulations','10000 simulations'},'location','northwest');
ylabel('number of tests not rejected'); xticklabels(fault_name); xtickangle(45); set(gca,'fontsize',11);
