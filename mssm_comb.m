%% Generate combined catalog from MSSM direct and MSSM adapt for use in PSHA codes

 %Combines MSSM and MSSM adapt into one catalog
 %Relative weightings based on cat_w as set in syncat_PSHA_MSSM_input
 %Can also make plot to compare MFD
clear all

load ('syncat_PSHA_MSSM_input','NumSimu','t_limit','cat_w','syncat_name','bg_mmin','max_bg',...
    'fault_width_weight','recurrence_type_weight');

addpath('syncat_MSSM');  addpath('syncat_adaptedMSSM'); addpath('syncat_bg');

load(char(syncat_name),'EQCAT'); load GR_TintiMulargia %load SSA-GEM catalog

%Event catalogs simplified to:
 % 1) event number
 % 2) occurrence time
 % 3) magnitude
 % 4) source_id
 % 5) width case (MSSM adapted only)
 % 6) char of G-R (MSSM adapted only)
 % 7) catalog

%either download MSSM_comb folder from Zenodo  
load MSSM_comb;

%or combine seperate catalogs by uncommenting below script
%{
MSSM_comb = [EQCAT(:,[1,2,5,7]),zeros(length(EQCAT),1),zeros(length(EQCAT),1),ones(length(EQCAT),1)];
    load ('MSSM_Catalog_Adapted_em','StochasticEventCatalog','dm');
    
    for i=1:4
        MSSM_comb=vertcat(MSSM_comb,[StochasticEventCatalog{i}(:,[1,2,4,6,5,7]),(i+1)*ones(length(StochasticEventCatalog{i}),1)]);
    end
    
NumSimu_fs=NumSimu*5; NumSimu_bg=NumSimu;

save('MSSM_comb','MSSM_comb')
%}
%% Freq mag plot comparison for all earthquake records

load syncat_bg %download from Zenodo

Mmin = bg_mmin; dm=0.005; NumSimu_fs=NumSimu*5; NumSimu_bg=NumSimu;

%Remove events outside Malawi and <7 in Malawi from bg_catalog
%(ie those actually used in PSHA)
rem_event=find((bg_catalog(:,4) >max_bg & bg_catalog(:,1)<0) | bg_catalog(:,1)>0);
bg_catalog(rem_event,:)=[];

%load data from MSSM fault-based syncat

all_mag_range_GR = ((Mmin-dm/2):0.05:8.2)';

allAnnualRate_comb=zeros(length(all_mag_range_GR),1);

allAnnualRate_fltw1_char=zeros(length(all_mag_range_GR),1); allAnnualRate_fltw1_exp=zeros(length(all_mag_range_GR),1); 
allAnnualRate_fltw2_exp=zeros(length(all_mag_range_GR),1); allAnnualRate_fltw2_char=zeros(length(all_mag_range_GR),1);
bg_AnnualRate=zeros(length(all_mag_range_GR),1); fault_syncatRate=zeros(length(all_mag_range_GR),1);

mo_indx1_exp=find(MSSM_comb(:,5) ==1 & MSSM_comb(:,6) ==2);
mo_indx1_char=find(MSSM_comb(:,5) ==1 & MSSM_comb(:,6) ==1);
mo_indx2_exp=find(MSSM_comb(:,5) ==2 & MSSM_comb(:,6) ==2);
mo_indx2_char=find(MSSM_comb(:,5) ==2 & MSSM_comb(:,6) ==1);

for ll = 1:length(all_mag_range_GR)
    
    bg_AnnualRate (ll,1) = length(find(bg_catalog(:,4) >= all_mag_range_GR(ll)))/(NumSimu_bg*t_limit); %Rate in Malawi from bg catalog
    
    allAnnualRate_comb(ll,1) = bg_AnnualRate (ll,1)+(length(find(MSSM_comb(:,3) >= all_mag_range_GR(ll)))/(NumSimu_fs*t_limit));%combined annual rate of events
    allAnnualRate_fltw1_exp(ll,1) = bg_AnnualRate (ll,1)+(length(find(MSSM_comb(mo_indx1_exp,3) >= all_mag_range_GR(ll)))/(NumSimu_fs*t_limit*cat_w(2)*fault_width_weight(1)*recurrence_type_weight(2)));%annual rate of char events, width=1
    allAnnualRate_fltw1_char(ll,1) = bg_AnnualRate (ll,1)+(length(find(MSSM_comb(mo_indx1_char,3) >= all_mag_range_GR(ll)))/(NumSimu_fs*t_limit*cat_w(2)*fault_width_weight(1)*recurrence_type_weight(1)));%annual rate of exp events, width=1
    allAnnualRate_fltw2_exp(ll,1) = bg_AnnualRate (ll,1)+(length(find(MSSM_comb(mo_indx2_exp,3) >= all_mag_range_GR(ll)))/(NumSimu_fs*t_limit*cat_w(2)*fault_width_weight(2)*recurrence_type_weight(2)));%annual rate of events where width=1
    allAnnualRate_fltw2_char(ll,1) = bg_AnnualRate (ll,1)+(length(find(MSSM_comb(mo_indx2_char,3) >= all_mag_range_GR(ll)))/(NumSimu_fs*t_limit*cat_w(2)*fault_width_weight(2)*recurrence_type_weight(1)));%annual rate of events where width=2

    fault_syncatRate (ll,1) = bg_AnnualRate (ll,1)+(length(find(EQCAT(:,5) >= all_mag_range_GR(ll)))/(NumSimu*t_limit));%Rate from MSSM direct catalog

end

%% Check Moment Rates

Mo_bg_catalog = 10.^(1.5*bg_catalog(:,4)+9.05);
MoRate_bg_catalog = sum(Mo_bg_catalog)/NumSimu_bg;

Mo_MSSM_Direct =  10.^(1.5*EQCAT(:,5)+9.05);
MoRate_MSSM_Direct = MoRate_bg_catalog +(sum(Mo_MSSM_Direct)/NumSimu);

Mo_StochasticEventCatalog_1_exp = 10.^(1.5*MSSM_comb(mo_indx1_exp,3)+9.05);
MoRate_StochasticEventCatalog_1_exp = MoRate_bg_catalog +(sum(Mo_StochasticEventCatalog_1_exp)/(NumSimu_fs*t_limit*cat_w(2)*fault_width_weight(1)*recurrence_type_weight(2)));

Mo_StochasticEventCatalog_1_char = 10.^(1.5*MSSM_comb(mo_indx1_char,3)+9.05);
MoRate_StochasticEventCatalog_1_char = MoRate_bg_catalog +(sum(Mo_StochasticEventCatalog_1_char)/(NumSimu_fs*cat_w(2)*t_limit*fault_width_weight(1)*recurrence_type_weight(1)));

Mo_StochasticEventCatalog_2_exp = 10.^(1.5*MSSM_comb(mo_indx2_exp,3)+9.05);
MoRate_StochasticEventCatalog_2_exp = MoRate_bg_catalog +(sum(Mo_StochasticEventCatalog_2_exp)/(NumSimu_fs*cat_w(2)*t_limit*fault_width_weight(2)*recurrence_type_weight(2)));

Mo_StochasticEventCatalog_2_char = 10.^(1.5*MSSM_comb(mo_indx2_char,3)+9.05);
MoRate_StochasticEventCatalog_2_char = MoRate_bg_catalog +(sum(Mo_StochasticEventCatalog_2_char)/(NumSimu_fs*cat_w(2)*t_limit*fault_width_weight(2)*recurrence_type_weight(1)));

Mo_all_MSSM =  10.^(1.5*MSSM_comb(:,3)+9.05);
MoRate_MSSM_all_MSSM = sum(Mo_all_MSSM)/NumSimu_fs;

MoRate_MSSM_comb =  MoRate_MSSM_all_MSSM  +MoRate_bg_catalog;

%% Make MFD Plot comparison

figure(2);

%Comparing different fault width scenarios

semilogy(all_mag_range_GR,allAnnualRate_fltw1_char,'r--v','LineWidth',1.5,'MarkerIndices',1:8:length(all_mag_range_GR),'MarkerSize',8); hold on; 
semilogy(all_mag_range_GR,allAnnualRate_fltw1_exp,'k--^','LineWidth',1.5,'MarkerIndices',1:8:length(all_mag_range_GR),'MarkerSize',8); hold on 
semilogy(all_mag_range_GR,allAnnualRate_fltw2_char,'r-','LineWidth',1.5); hold on; 
semilogy(all_mag_range_GR,allAnnualRate_fltw2_exp,'k-','LineWidth',1.5); hold on 
semilogy(all_mag_range_GR,fault_syncatRate,'--*','LineWidth',1.5,'MarkerIndices',4:8:length(all_mag_range_GR),'MarkerSize',8); hold on;
semilogy(all_mag_range_GR,allAnnualRate_comb,'b-','LineWidth',2); hold on; %MFD from event catalog
  axis([4.5 8.1 10^-5 5*10^0]);
legend({['Length-limited {\it W}, Char,' newline 'M{_0} : ',num2str(MoRate_StochasticEventCatalog_1_char,2),' Nm/yr'],...
    ['Length-limited {\it W}, G-R,' newline 'M{_0}: ',num2str(MoRate_StochasticEventCatalog_1_exp,2),' Nm/yr'],...
    ['Layer-limited {\it W}, Char,' newline 'M{_0}: ', num2str(MoRate_StochasticEventCatalog_2_char,2),' Nm/yr'],...
    ['Layer-limited {\it W}, G-R,' newline 'M{_0}: ' , num2str(MoRate_StochasticEventCatalog_2_exp,2),' Nm/yr'],...
    ['MSSM Direct,' newline 'M{_0}: ',num2str(MoRate_MSSM_Direct,2),' Nm/yr'],...
    ['Combined,' newline 'M{_0}: ',num2str(MoRate_MSSM_comb,2),' Nm/yr']},...
    'Location','southwest','FontSize',10);hold on;
set(gca,'fontsize',13);
xlabel('Magnitude'); ylabel('Annual frequency of exceedance'); grid on; axis square;



%% Plot moment rate tests for MSSM Direct and MSSM adapted comparison

load syncat_comparison; load MSSM_adapted_comparison

figure(104);

subplot(1,2,1)

plot([13.5 17.5],[13.5 17.5],'k--','LineWidth',1.5); hold on 
plot(log10(num_sec1(:,26)),log10(syncat_sec_mr),'rx','MarkerSize',12,'LineWidth',1.2); hold on
plot(log10(num_fault(:,25)),log10(syncat_fault_mr),'bx','MarkerSize',12,'LineWidth',1.2); hold on
plot(log10(num_multi_fault(:,15)),log10(syncat_multi_fault_mr),'x','MarkerSize',12,'Color',[0 0.7 0],'LineWidth',1.2);

t1 = title('(a)','fontweight','normal');
set(t1, 'horizontalAlignment', 'left','units','normalized');
h1 = get(t1, 'position');
set(t1, 'position', [-0.2 h1(2) h1(3)]); set(gca,'FontSize',13)

axis([13.5 17.5 13.5 17.5]); axis square; hold on;
xlabel('MSSM Moment Rate (log Nm/yr)'),ylabel(['Stochastic Event Catalog' newline 'Moment Rate (log Nm/yr)']); hold on;
legend('y=x',['Section Sources n = ' num2str(length(num_sec1))],['Fault Sources n=' num2str(length(num_fault))],...
    ['Multi-fault Sources n =' num2str(length(num_multi_fault))] ,'Location','southeast');
grid on; set(gca,'fontsize',13)


subplot(1,2,2)
plot([14.5 17.5],[14.5 17.5],'k--'); hold on
%plots for case where ss=1 and ss=3 (i.e char cases with different widths)
plot(log10(num_MSSM(:,13)),log10(MoRate_StochasticEventCatalog{1}(:,1)),'rx',...
log10(num_MSSM(:,14)),log10(MoRate_StochasticEventCatalog{4}(:,1)),'kx','MarkerSize',12,'LineWidth',1.5);

t1 = title('(b)','fontweight','normal');
set(t1, 'horizontalAlignment', 'left','units','normalized');
h1 = get(t1, 'position');
set(t1, 'position', [-0.2 h1(2) h1(3)]); set(gca,'FontSize',13)


axis([14.5 17.5 14.5 17.5]); axis square; hold on;
xlabel('MSSM Moment Rate (log Nm/yr)'),ylabel(['Stochastic Event Catalog' newline 'Moment Rate (log Nm/yr)']); hold on;
legend('y=x','Length-limited {\it W}, Char','Layer-limited {\it W}, G-R','Location','southeast');
grid on; set(gca,'fontsize',13)

set(gcf,'Position',[681 566 821 413]);

%% Sample and compare Moment Rate in 50 year samples of MSSM_comb

%randomly assign simulation cycle to bg_catalog
bg_catolog_simu= NumSimu.*rand(length(bg_catalog),1); 
simu_count = 1; mo_test = {}; clock_old = cputime;

for ii=1:5
    
 %Sample each catalog dy duration of assessed GEM catalog (50 years), taking only simulation cycle and mag
 tmp_catalog=MSSM_comb(find(MSSM_comb(:,7)==ii),2:3);
 tmp_catalog=vertcat(tmp_catalog,[bg_catolog_simu bg_catalog(:,4)]); 
 
 mo_test{ii}=zeros(NumSimu/TotalYear,3); 
 
 for jj=1:NumSimu/TotalYear
       
   if jj == 10000*simu_count
        clock_new = cputime;
        disp(['Simulation cycle is: ',num2str(jj),' & Catalog = ',num2str(ii),' & Required time is: ',num2str(clock_new-clock_old)]);
        clock_old = clock_new;
        simu_count = simu_count + 1;
   end  
     
     jj_indx=jj*TotalYear;
     indx=find(tmp_catalog(:,1)<jj_indx & tmp_catalog(:,1)>(jj_indx-TotalYear));
     mo_test{ii}(jj,1)=sum(10.^(1.5.*tmp_catalog(indx,2)+9.05))/TotalYear;%find moment rate of each sample
     mo_test{ii}(jj,2)=max(tmp_catalog(indx,2));%find Mmax for each sample
     mo_test{ii}(jj,3)=length(find(tmp_catalog(indx,2)>=6.0));%find number of events > Mw6 forof each sample
 end
  
 simu_count = 1;
 
end

save('mo_test','mo_test');

%% Plot analysis

load geo_mo_rate 
load mo_test; %uncomment if mo_test already run

mo_test_mat = cell2mat(mo_test');

%Find observed moment rate, mmax and events >6 for 1965-2015 from GEM catalog 
indx=find(EAcatalog_R1(:,10)>bg_mmin);
ins_mo_test = [sum(10.^(1.5.*EAcatalog_R1(indx,10)+9.05))/TotalYear, max(EAcatalog_R1(indx,10)), ...
         length(find(EAcatalog_R1(:,10)>=6.0))];
     
for jj=1:length(mo_test)
    
    % Percent change that moment rate in MSSM_comb is lower than observed
    mo_rate_prob(jj)=length(find(mo_test{jj}(:,1)<=ins_mo_test(:,1)))/length(mo_test{jj}(:,1))*100;
    % Percent change that mmax in MSSM_comb is lower than observed
    m_max_prob(jj)=length(find(mo_test{jj}(:,2)<=ins_mo_test(:,2)))/length(mo_test{jj}(:,1))*100;
    
end

total_mo_rate_prob=length(find(mo_test_mat(:,1)<=ins_mo_test(:,1)))/length(mo_test_mat(:,1))*100;
total_mmax_prob=length(find(mo_test_mat(:,2)<=ins_mo_test(:,2)))/length(mo_test_mat(:,2))*100;    
  
%FIGURE OPTION 1: PLOTTING MOMENT DISTRIBUTION WITH INSTRUMENTAL RECORD

subplot(1,2,1)

marker_I = ['v','*','v','*'];
line_S = ["-","-."];

for i=1:length(mo_test)
    pd_i = fitdist(log10(mo_test{i}(:,1)),'kernel');
    x_i = log10(min(mo_test{i}(:,1))):0.05:log10(max(mo_test{i}(:,1)));
    y_i = pdf(pd_i,x_i); ylim=1.5; x_inc =0.05;
    if i==1
    plot(x_i,y_i,'Marker','o','color',[0.4 0.4 0.4],'linewidth',1.3,'MarkerIndices',[1:round(length(x_i)/16):length(x_i)-1],'MarkerFaceColor',[1 1 1]);hold on
    else
    plot(x_i,y_i,line_S(floor(i/2)),'linewidth',1.3,'color',[0 0 0],'Marker',marker_I(i-1),...
        'MarkerIndices',[i:round(length(x_i)/16):length(x_i)-1],'MarkerFaceColor',[1 1 1]);hold on
    end
end
        
plot_data2=log10(mo_test_mat(:,1));
plot_data3=log10(ins_mo_test(:,1)); xlabel_text='log10 Moment Rate';   
prob_text = horzcat(mo_rate_prob,total_mmax_prob); title_text = ('(a)');
observed_text1 = ['Observed seismic' newline 'moment rate'];
observed_text2 = 'GSRM v2.1';
observed_text3 = 'SSA-GSRM v1.0';
             
hold on

pd = fitdist(plot_data2,'kernel');
x = min(plot_data2):x_inc:max(plot_data2);
y = pdf(pd,x);
plot(x,y,'b-','linewidth',1.9);hold on

%Plot observed moment rate

plot([plot_data3 plot_data3],[0 ylim],'r-.','linewidth',2);

%Plot moment rates from geodesy
plot([log10(mo_rate_opt(:,1)) log10(mo_rate_opt(:,1))],[0 1.5],'color',[0 0.5 0],'linestyle','-.','linewidth',2);
plot([log10(mo_rate_opt(:,2)) log10(mo_rate_opt(:,2))],[0 1.5],'m-.','linewidth',2);

xlabel(xlabel_text); ylabel('f(x)'); set(gca,'fontsize',13); axis square

t1 = title('(a)','fontweight','normal');
set(t1, 'horizontalAlignment', 'left','units','normalized');
h1 = get(t1, 'position');
set(t1, 'position', [-0.2 h1(2) h1(3)]); 

legend(['Direct MSSM: ',num2str(prob_text(:,1),2),'%'],['Length-limited {\it W}' newline 'Char: ' num2str(prob_text(:,2),2) '%'],...
        ['Length-limited {\it W},' newline 'G-R: ' num2str(prob_text(:,3),1) '%'],['Layer-limited {\it W},' newline 'Char: ' num2str(prob_text(:,4),2) '%'],...
        ['Layer-limited {\it W},' newline 'G-R: ' num2str(prob_text(:,5),1) '%'],['Combined ' num2str(prob_text(:,6),2) '%'], observed_text1,observed_text2,observed_text3,'fontsize',9);

%Plot instrumental MFD
subplot(1,2,2)

semilogy(all_mag_range_GR,allAnnualRate_comb,'b-','LineWidth',2); hold on; %MFD from event catalog
semilogy(mag_check,AnnualCumNum,'r*--','LineWidth',1.5); hold on;
semilogy(4:0.1:8.1,10.^(Stat_log10N(:,4)),'k-',4:0.1:8.1,10.^(Stat_log10N(:,[3 5])),'k--','LineWidth',1);
        axis square; xlabel('Magnitude'); ylabel('Annual frequency of exceedance'); 
legend({'Combined','SSA-GEM Catalog 1965-2015','Fitted G-R relationship','16-th/84-th percentile'},...
    'Location','southwest','FontSize',9.5);
set(gca,'fontsize',13); grid on;       
axis([4.5 8.1 10^-5 5*10^0]); axis square; hold on;

t2=title('(b)','fontsize',14,'fontweight','normal');
set(t2, 'horizontalAlignment', 'left');
set(t2, 'units', 'normalized');
h1 = get(t2, 'position');
set(t2, 'position', [-0.2 h1(2) h1(3)]);

set(gcf,'Position',[440 314 1122 483]);

%% FIGURE OPTION 2
%Fit Kernel density function to sampled Moment Rate
%Change kk to 2 for plotting Mmax or 3 for number of events Mw >6

figure(999);

for kk=1:2
    
    subplot(1,2,kk)

    if kk==1
        
        for i=1:length(mo_test)
            pd_i = fitdist(log10(mo_test{i}(:,kk)),'kernel');
            x_i = log10(min(mo_test{i}(:,kk))):0.05:log10(max(mo_test{i}(:,kk)));
            y_i = pdf(pd_i,x_i); ylim=1.5; x_inc =0.05;
            plot(x_i,y_i,'linewidth',1.3);hold on
        end
        
        plot_data2=log10(mo_test_mat(:,kk));
        plot_data3=log10(ins_mo_test(:,kk)); xlabel_text='log10 Moment Rate';   
        prob_text = horzcat(mo_rate_prob,total_mmax_prob); title_text = ('(a)');
        observed_text1 = ['Observed seismic' newline 'moment rate'];
        observed_text2 = 'GSRM v2.1';
        observed_text3 = 'SSA-GSRM v1.0';
        
    else
        
        for i=1:length(mo_test)
            pd_i = fitdist(mo_test{i}(:,kk),'kernel');
            x_i = min(mo_test{i}(:,kk)):0.05:max(mo_test{i}(:,kk));
            y_i = pdf(pd_i,x_i); ylim=1.5; x_inc =0.05;
            plot(x_i,y_i,'linewidth',1.3);hold on
        end
        
        plot_data2=mo_test_mat(:,kk); observed_text = 'Observed M_{Max}';
        plot_data3=ins_mo_test(:,kk); xlabel_text='M_{Max}';
        prob_text = horzcat(m_max_prob,total_mo_rate_prob); title_text = ('(b)');
    end 
        
hold on

pd = fitdist(plot_data2,'kernel');
x = min(plot_data2):x_inc:max(plot_data2);
y = pdf(pd,x);
plot(x,y,'k-','linewidth',1.5);hold on

%Plot observed moment rate

plot([plot_data3 plot_data3],[0 ylim],'r--','linewidth',2);

%Plot moment rates from geodesy
if kk==1
    plot([log10(mo_rate_opt(:,1)) log10(mo_rate_opt(:,1))],[0 1.5],'b--','linewidth',2);
    plot([log10(mo_rate_opt(:,2)) log10(mo_rate_opt(:,2))],[0 1.5],'m--','linewidth',2);
end

xlabel(xlabel_text); ylabel('f(x)'); set(gca,'fontsize',13); axis square


    t1 = title(title_text,'fontweight','normal');
    set(t1, 'horizontalAlignment', 'left','units','normalized');
    h1 = get(t1, 'position');
    set(t1, 'position', [-0.2 h1(2) h1(3)]); 


    if kk==1
        legend(['Direct MSSM: ',num2str(prob_text(:,1),2),'%'],['length-limited {\it W}' newline 'Char: ' num2str(prob_text(:,2),2) '%'],...
        ['Length-limited {\it W},' newline 'G-R: ' num2str(prob_text(:,3),1) '%'],['Layer-limited {\it W},' newline 'Char: ' num2str(prob_text(:,4),2) '%'],...
        ['Layer-limited {\it W},' newline 'G-R: ' num2str(prob_text(:,5),1) '%'],['Combined ' num2str(prob_text(:,6),2) '%'], observed_text1,observed_text2,observed_text3,'fontsize',9);

    else
    
        legend(['Direct MSSM: ',num2str(prob_text(:,1),2),'%'],['length-limited {\it W}' newline 'Char: ' num2str(prob_text(:,2),2) '%'],...
        ['Length-limited {\it W},' newline 'G-R: ' num2str(prob_text(:,3),1) '%'],['Layer-limited {\it W},' newline 'Char: ' num2str(prob_text(:,4),2) '%'],...
        ['Layer-limited {\it W},' newline 'G-R: ' num2str(prob_text(:,5),1) '%'],['Combined ' num2str(prob_text(:,6),2) '%'], observed_text,'fontsize',9);
    end
end

set(gcf,'Position',[440 412 907 385]);





