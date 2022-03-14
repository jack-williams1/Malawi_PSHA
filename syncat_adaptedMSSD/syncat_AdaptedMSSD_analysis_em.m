%% PERFORM ANALYSIS OF STOCHASTIC EVENT CATALOG OF MSSD_ADAPTED CATALOG

load MSSD_Catalog_Adapted_em

% StochasticEventCatalog
  % 1    ) Overall cycle
  % 2    ) Simulation cycle
  % 3    ) Occurrence time
  % 4    ) Earthquake magnitude
  % 5    ) Fault model (1 = partial width versus 2 = full width)
  % 6    ) MSSD ID
  % 7    ) Recurrence type (1 = characteristic versus 2 = exponential)
  % 8    ) Reccurrece parameter

mydir  = pwd; idcs   = strfind(mydir,'/');
addpath(mydir(1:idcs(end)-1));
addpath([mydir(1:idcs(end)-1) '/syncat_bg']); addpath([mydir(1:idcs(end)-1) '/syncat_MSSD']);
addpath([mydir(1:idcs(end)-1) '/misc_functions']);

load MSSD_sources; load ('EQCAT_test','EQCAT','num_sec1');
load ('syncat_PSHA_MSSD_input','weight_SWM');

%% Result check

%syncat opt
%1) fault_width=1 & char
%2) fault_width=1 & GR
%3) fault_width=2 & char
%4) fault_width=2 & GR

% Seismic moment for whole catalog
for ss=1:length(syncat_opt)

Mo_StochasticEventCatalog_all = 10.^(1.5*StochasticEventCatalog{ss}(:,4)+9.05);
MoRate_StochasticEventCatalog_all(ss) = sum(Mo_StochasticEventCatalog_all)/(NumSimu*t_limit);
end

%% Seismic moment for each fault and for different width cases:

f_indx={}; f_indx1={}; f_indx2={};
Mo_StochasticEventCatalog={};


for ss=1:length(syncat_opt)
    
    fault_width=syncat_opt(ss,1);
    MoRate_StochasticEventCatalog{ss}=zeros(length(num_MSSD),2);
    for ff=1:length(num_MSSD)
    
    f_indx{ff}=find(StochasticEventCatalog{ss}(:,6) == num_MSSD(ff,1)); %Index all events for fault
    Mo_StochasticEventCatalog = 10.^(1.5*StochasticEventCatalog{ss}(f_indx{ff},4)+9.05);
    MoRate_StochasticEventCatalog{ss}(ff,1) = sum(Mo_StochasticEventCatalog)/(NumSimu*t_limit);

    %Theoritical moment rates
    if syncat_opt(ss,2)==1
        MoRate_StochasticEventCatalog{ss}(ff,2)=YC85SeismicMoment_char{ff}(5,fault_width);
    else
        MoRate_StochasticEventCatalog{ss}(ff,2)=YC85SeismicMoment_exp{ff}(5,fault_width);
    end %end if statement
    end %end ff loop
    clear fault_width f_indx  Mo_StochasticEventCatalog 
end

%% Moment rate comparison plots %Saves variable MSSD_adapted_comparison to be used in MSSD_comb

figure(1);
%Compare to MSSD moment rate estimates
subplot(1,2,1)
plot([14.5 17.5],[14.5 17.5],'k--'); hold on
%plots for case where ss=1 and ss=3 (i.e char cases with different widths)
plot(log10(num_MSSD(:,13)),log10(MoRate_StochasticEventCatalog{1}(:,1)),'bx',...
log10(num_MSSD(:,14)),log10(MoRate_StochasticEventCatalog{3}(:,1)),'rx','MarkerSize',12,'LineWidth',1.5);

axis([14.5 17.5 14.5 17.5]); axis square; hold on;
xlabel('MSSD Moment Rate (log Nm/yr)'),ylabel('Stochastic Event Catalog Moment Rate (log Nm/yr)'); hold on;
legend('y=x','length-limited width','layer-limited width','Location','southeast');
grid on; set(gca,'fontsize',13)

save('MSSD_adapted_comparison','MoRate_StochasticEventCatalog','num_MSSD');

subplot(1,2,2)
%Compare to Y&C85 moment rate estimates when fault width limited by fault
%length
plot([14.5 17.5],[14.5 17.5],'k--'); hold on

plot(log10(MoRate_StochasticEventCatalog{1}(:,2)),log10(MoRate_StochasticEventCatalog{1}(:,1)),'bx',...
log10(MoRate_StochasticEventCatalog{2}(:,2)),log10(MoRate_StochasticEventCatalog{2}(:,1)),'kx','MarkerSize',12,'LineWidth',1.5);

axis([14.5 17.5 14.5 17.5]); axis square; hold on;
xlabel('Y&C 85 Reccurrence Model Moment Rate (log Nm/yr)'),ylabel('Stochastic Event Catalog Moment Rate (log Nm/yr)'); hold on;
legend('y=x','char distribution','G-R distribution','Location','southeast');
grid on; set(gca,'fontsize',13)

set(gcf,'Position', [681 533 904 446]);

%% Freq mag plot comparison for all earthquake records

load syncat_bg

%Remove events outside Malawi 
rem_event=find(bg_catalog(:,1)<0);

%if only want to consider events <7 from bg_catalog
%(ie those actually used in PSHA)
%max_bg = xlsread(inputexcel,1,'D43');
%rem_event=find(bg_catalog(:,4) >max_bg & bg_catalog(:,1)<0);

bg_catalog(rem_event,:)=[]; Mmin=4.5;

%load data from MSSD fault-based syncat

all_mag_range_GR = ((Mmin-dm/2):0.05:8.5)'; bg_AnnualRate = zeros(length(all_mag_range_GR),1);

MSSD_adaptedAnnualRate=zeros(length(all_mag_range_GR),length(syncat_opt)); 


for ll = 1:length(all_mag_range_GR)
    bg_AnnualRate (ll,1) = length(find(bg_catalog(:,4) >= all_mag_range_GR(ll)))/(NumSimu*t_limit); %Rate in Malawi from bg catalog
    for ss=1:length(syncat_opt)
    MSSD_adaptedAnnualRate(ll,ss) = length(find(StochasticEventCatalog{ss}(:,4) >= all_mag_range_GR(ll)))/(NumSimu*t_limit);%combined annual rate of events
    MSSD_adaptedAnnualRate_comb(ll,ss) = MSSD_adaptedAnnualRate(ll,ss)+bg_AnnualRate(ll);
    end
    
    fault_syncatRate (ll,1) = length(find(EQCAT(:,5) >= all_mag_range_GR(ll)))/(NumSimu*t_limit);%Rate from MSSD direct catalog

end

%% Make MFD Plot comparison

figure(2);

%Comparing all different fault-based catalogs
semilogy(all_mag_range_GR,MSSD_adaptedAnnualRate(:,1),'k-','LineWidth',1.2); hold on; %MFD from event catalog
semilogy(all_mag_range_GR,MSSD_adaptedAnnualRate(:,2),'k--x','LineWidth',1.2); hold on; %MFD from event catalog
semilogy(all_mag_range_GR,MSSD_adaptedAnnualRate(:,3),'-','Color',[0.6 0.6 0.6],'LineWidth',1.2); hold on; %MFD from event catalog
semilogy(all_mag_range_GR,MSSD_adaptedAnnualRate(:,4),'--x','Color',[0.6 0.6 0.6],'LineWidth',1.2); hold on; %MFD from event catalog
semilogy(all_mag_range_GR,fault_syncatRate,'b-','LineWidth',1.2); hold on; %MFD from event catalog
semilogy(all_mag_range_GR,bg_AnnualRate,'r-','LineWidth',1.2); hold on; %MFD from event catalog
set(gca,'fontsize',13);  axis([5.5 8.5 10^-5 10^-1]); legend('length-limited w, char','length-limited w, G-R',...
    'layer-limited w, char','layer-limited w, GR','MSSD Direct','background','Location','southwest');hold on;
xlabel('Magnitude'); ylabel('Annual frequency'); grid on; axis square;


%% Freq mag plot for particular fault or multifault

ff=314; %select fault or multi_fault based on ID
indx=find(num_MSSD(:,1)==ff);
flt_ruptures={}; sec_ruptures={};


if ff>600
    %finf multifault ruptures in EQCAT
    m_flt_id = find(ff==num_multi_fault(:,1));
    m_flt_ruptures=find(EQCAT(:,7)==ff);
    
    %find participating fault of multifault system
    flt_id=num_fault(find(ff==num_fault(:,23)),1);
    
    %find fault ruptures in MSSD
    for i=1:length(flt_id)
        flt_indx(i)=find(num_fault(:,1)==flt_id(i));
        flt_ruptures{i}=find(EQCAT(:,7)==flt_id(i));
        if num_fault(flt_indx(i),13)>0%if fault is segmented
        weighting_f(i) = 1/3;%weighting between multi_fault and fault,ruptures
        else
        weighting_f(i) = 0.5;%weighting between multi_fault and fault,ruptures
        end
    end
    
    %find sec in MSSD
    sec_id=num_sec1(find(ff==num_sec1(:,18)),1);  
    if isempty(sec_id)==0;
        for i=1:length(sec_id)
            %Index num_sec1 where no linking sections
            sec_indx(i)=find(num_sec1(:,1)==sec_id(i));
            sec_ruptures{i}=find(EQCAT(:,7)==sec_id(i));
        end 
        weighting_s = 1/3;%weighting between multi_fault ,fault, and sec ruptures
        weighting = 1/3;%weighting between multi_fault ,fault, and sec ruptures
    else
        weighting = 0.5;%weighting between multi_fault ,fault, and sec ruptures
    end
    
else
    
    flt_ruptures{1}=find(EQCAT(:,7)==ff);
    %find if participating sections in fault system
    flt_indx=find(ff==num_fault(:,1));
    if num_fault(flt_indx,13)>0
        sec_id=num_sec1(find(ff==num_sec1(:,15)),1); 
        for i=1:length(sec_id)
        %Index num_sec1 where no linking sections
        sec_indx(i)=find(num_sec1(:,1)==sec_id(i));
        sec_ruptures{i}=find(EQCAT(:,7)==sec_id(i));
        end
        weighting_s = 0.5;%weighting between fault and sec ruptures
        weighting = 0.5;
    else
        weighting =1; sec_id =[];
    end
end

sourceRate={};

%Index events from ff in StochasticEventCatalog
for ss=1:length(StochasticEventCatalog)
source_events=find(StochasticEventCatalog{ss}(:,6) == ff); 
sourceRate{ss}=zeros(length(mag_range_GR{indx}),1);
%Index events based on recurrence model and fault width
    for ii = 1:length(mag_range_GR{indx})
       sourceRate{ss}(ii,1) = length(find(StochasticEventCatalog{ss}(source_events,4) >= mag_range_GR{indx}(ii)))/(NumSimu*t_limit);%annual rate of events where width=1
    end
end

%Find model recurrence models for source ff
weighted_GR_char1 = zeros(length(mag_range_GR{indx}),1);
weighted_GR_exp1  = zeros(length(mag_range_GR{indx}),1);  

weighted_GR_char2 = zeros(length(mag_range_GR{indx}),1);
weighted_GR_exp2  = zeros(length(mag_range_GR{indx}),1); 

for jj = 1:height(RecurrenceVar{indx})
   weighted_GR_char1 = weighted_GR_char1 + RecurrenceVar{indx}(jj,5)*fm_char_var1{indx}(:,3*jj);
   weighted_GR_exp1  = weighted_GR_exp1  + RecurrenceVar{indx}(jj,5)*fm_exp_var1{indx}(:,3*jj);
   
   weighted_GR_char2 = weighted_GR_char2 + RecurrenceVar{indx}(jj,5)*fm_char_var2{indx}(:,3*jj);
   weighted_GR_exp2  = weighted_GR_exp2  + RecurrenceVar{indx}(jj,5)*fm_exp_var2{indx}(:,3*jj);
end

% Mfd comparison between Stochastic Event Catalog and Recurrence Models

figure(4);
%Plot for intermediate estimates
%plot width =1
semilogy(mag_range_GR{indx},fm_char_var1{indx}(:,((length(RecurrenceVar{indx})/2)+0.5)*3),'r-',mag_range_GR{indx},fm_exp_var1{indx}(:,((length(RecurrenceVar{indx})/2)+0.5)*3),'b-' ,'LineWidth',1.5); hold on; 
%plot width =2
semilogy(mag_range_GR{indx},fm_char_var2{indx}(:,((length(RecurrenceVar{indx})/2)+0.5)*3),'m-',mag_range_GR{indx},fm_exp_var2{indx}(:,((length(RecurrenceVar{indx})/2)+0.5)*3),'k-','LineWidth',1.5); hold on;
%MFD from event catalog
for ss=1:length(StochasticEventCatalog)
    semilogy(mag_range_GR{indx},sourceRate{ss},'Color',[0.3 0.3 0.3],'LineStyle','--','LineWidth',1.8);  
    
end
      
legend(['Char events:' newline 'length limited width'],['G-R events:' newline 'length limited width'],...
       ['Char events:' newline 'layer limited width'],['G-R events:' newline 'layer limited width'],['Simulated Events'], 'Location','southwest'); hold on;
set(gca,'fontsize',13); axis([5 num_MSSD(indx,12)+0.2 10^-6 10^-2]); hold on;
xlabel('Magnitude'); ylabel('Annual frequency'); grid on; axis square;

figure(5);
%Plot for weighted average from recurrence models
% Note weighted average cuts out at mmax-mshift (i,e lower bound of Mmax)
%plot width =1
semilogy(mag_range_GR{indx},weighted_GR_char1,'r-',mag_range_GR{indx},weighted_GR_exp1,'b-' ,'LineWidth',1.5); hold on; 
%plot width =2
semilogy(mag_range_GR{indx},weighted_GR_char2,'m-',mag_range_GR{indx},weighted_GR_exp2,'k-','LineWidth',1.5); hold on;
%MFD from event catalog
for ss=1:length(StochasticEventCatalog)
    semilogy(mag_range_GR{indx},sourceRate{ss},'Color',[0.3 0.3 0.3],'LineStyle','--','LineWidth',1.8);  
    
end      
legend(['Char events:' newline 'length limited width'],['G-R events:' newline 'length limited width'],...
       ['Char events:' newline 'layer limited width'],['G-R events:' newline 'layer limited width'], 'Location','southwest'); hold on;
set(gca,'fontsize',13); axis([5 num_MSSD(indx,12)+0.2 10^-6 10^-2]); hold on;
xlabel('Magnitude'); ylabel('Annual frequency'); grid on; axis square;

%% Plot adapted_MSSD results with MSSD itself and direct MSSD eqcat catalog 
%Only considers fault width = 1 (for concistency with MSSD database)

figure(6);
a1=semilogy(mag_range_GR{indx},fm_char_var1{indx}(:,((length(RecurrenceVar{indx})/2)+0.5)*3),'r--',mag_range_GR{indx},fm_exp_var1{indx}(:,((length(RecurrenceVar{indx})/2)+0.5)*3),'k--' ,'LineWidth',1.5); hold on; 
%plot width =2
a2=semilogy(mag_range_GR{indx},fm_char_var2{indx}(:,((length(RecurrenceVar{indx})/2)+0.5)*3),'r-',mag_range_GR{indx},fm_exp_var2{indx}(:,((length(RecurrenceVar{indx})/2)+0.5)*3),'k-','LineWidth',1.5); hold on;

if ff>600 %For multi_fault, plot results from MSSD and eqcat

    m1=semilogy(num_multi_fault(m_flt_id,10),(1/num_multi_fault(m_flt_id,13)),'k*'); hold on; %Plot M-F from MSSD for multi_fault
    
    for i=1:length(flt_id)
        f1(i)=semilogy(num_fault(flt_indx(i),18),(1/num_fault(flt_indx(i),21)),'LineWidth',1.3,'*','Markersize',8,'MarkerEdgeColor',[0 0.5 0]); hold on; %Plot M-F from MSSD for fault
      if i>1 %avoid duplication for legend indexing
         set(get(get(p1(i),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
      end
    end
        
    if isempty(sec_id)==0 %i.e. fault is segmented
        
    for i=1:length(sec_id)
        s1(i)=semilogy(num_sec1(sec_indx(i),21),(1/num_sec1(sec_indx(i),24)),'LineWidth',1.3,'*','Markersize',8,'MarkerEdgeColor',[0.4 0.4 0.4]); hold on; %Plot M-F from MSSD for sec
      if i>1  %avoid duplication for legend indexing
         set(get(get(s1(i),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
      end
        
    end
        legend([s1 f1 m1 a1(1) a1(2) a2(1) a2(2)],{['MSSD MultiFault', 'MSSD Fault','MSSD Section','Adapted char events:' newline 'length-limited width'],['Adapted G-R events:' newline 'length-limited width'],...
            ['Adapted char events:' newline 'layer-limited width'],['adapted G-R events:' newline 'layer-limited width']},...
            'Location','northeast','fontsize',12); hold on;
    else
         legend([m1 f1 a1(1) a1(2) a2(1) a2(2)],{['MSSD MultiFault'], ['MSSD Fault'],['Adapted char events:' newline 'length-limited width'],['Adapted G-R events:' newline 'length-limited width'],...
             ['Adapted char events:' newline 'later-limited width'],['Adapted G-R events:' newline 'layer-limited width']},...
            'Location','northeast','fontsize',12); hold on;
    end
    
else %fault is not part of multi_fault system
    
    f1=semilogy(num_fault(flt_indx,18),(1/num_fault(flt_indx,21)),'b*','LineWidth',1.3,'Markersize',8); hold on; %Plot M-F from MSSD for fault
    
    if isempty(sec_id)==0 %i.e. fault is segmented
        
    for i=1:length(sec_id)
        s1=semilogy(num_sec1(sec_indx(i),21),(1/num_sec1(sec_indx(i),24)),'*','LineWidth',1.3,'Markersize',8,'MarkerEdgeColor',[0.4 0.4 0.4]); hold on; %Plot M-F from MSSD for sec
        
        legend([s1 f1 a1(1) a1(2) a2(1) a2(2)],{['MSSD Section'],['MSSD Fault'],['Adapted char events:' newline 'length-limited width'],['Adapted G-R events:' newline 'length-limited width'],...
            ['Adapted char events:' newline 'layer-limited width'],['Adapted G-R events:' newline 'layer-limited width']},...
             'Location','northeast','fontsize',11); hold on;
    end
    else
        legend([f1 a1(1) a1(2) a2(1) a2(2)],{['MSSD Fault'],['Adapted char events:' newline 'length-limited width'],['Adapted G-R events:' newline 'length-limited width'],...
            ['Adapted char events:' newline 'layer-limited width'],['Adapted G-R events:' newline 'layer-limited width']},...
            'Location','northeast','fontsize',11); hold on;
    end 

end

axis([5 num_MSSD(indx,12)+0.2 10^-6 10^-2]); xlabel('Magnitude'); ylabel('Annual occurrence rate'); axis square; grid on; set(gca,'fontsize',15)

%% Plot all theortical model curves for the MSSD source
figure(7)

subplot(121); 

p1=semilogy(mag_range_GR{indx},fm_char_var1{indx}(:,3:3:end),'Color',[0.5 0.5 0.5],'LineWidth',1); hold on
p2=semilogy(mag_range_GR{indx},weighted_GR_char1,'k-','LineWidth',1.5); hold on; 

xlabel('Magnitude'); ylabel('Annual frequency of exceedance');
legend([p1(1) p2],'Individual','Weigted Average','Location','Southwest');
axis([6 8 10^-6 10^-3]); axis square; grid on
t1 = title('(a)');
set(t1, 'horizontalAlignment', 'left','units','normalized');
h1 = get(t1, 'position');
set(t1, 'position', [-0.2 h1(2) h1(3)]); set(gca,'FontSize',13)

subplot(122);

p1=semilogy(mag_range_GR{indx},fm_exp_var1{indx}(:,3:3:end),'Color',[0.5 0.5 0.5],'LineWidth',1); hold on
p2=semilogy(mag_range_GR{indx},weighted_GR_exp1,'k-','LineWidth',1.5); hold on; 

xlabel('Magnitude'); ylabel('Annual frequency of exceedance');
legend([p1(1) p2],'Individual','Weigted Average','Location','Southwest');
axis([6 8 10^-6 10^-3]); axis square; grid on
t1 = title('(b)');
set(t1, 'horizontalAlignment', 'left','units','normalized');
h1 = get(t1, 'position');
set(t1, 'position', [-0.2 h1(2) h1(3)]); set(gca,'FontSize',13)

set(gcf,'Position', [681 629 726 350]);