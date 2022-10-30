%% PERFORM ANALYSIS OF STOCHASTIC EVENT CATALOG OF MSSM_ADAPTED CATALOG

load MSSM_Catalog_Adapted_em

% StochasticEventCatalog
  % 1    ) Overall cycle
  % 2    ) Simulation cycle
  % 3    ) Occurrence time
  % 4    ) Earthquake magnitude
  % 5    ) Fault model (1 = partial width versus 2 = full width)
  % 6    ) MSSM ID
  % 7    ) Recurrence type (1 = characteristic versus 2 = exponential)
  % 8    ) Reccurrece parameter

mydir  = pwd; idcs   = strfind(mydir,'/');
addpath(mydir(1:idcs(end)-1));
addpath([mydir(1:idcs(end)-1) '/syncat_bg']); addpath([mydir(1:idcs(end)-1) '/syncat_MSSM']);
addpath([mydir(1:idcs(end)-1) '/misc_functions']);

load MSSM_sources; load ('EQCAT_test','EQCAT','num_sec1');
load ('syncat_PSHA_MSSM_input','weight_SWM','t_limit','NumSimu');
MSSM_adapt=readtable([mydir(1:idcs(end)-1) '/MSSM_sourcecalcs/MSSM.xlsx'],'sheet','MSSM_AdaptedSources');

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
    MoRate_StochasticEventCatalog{ss}=zeros(length(num_MSSM),2);
    for ff=1:length(num_MSSM)
    
    f_indx{ff}=find(StochasticEventCatalog{ss}(:,6) == num_MSSM(ff,1)); %Index all events for fault
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

%% Moment rate comparison plots %Saves variable MSSM_adapted_comparison to be used in MSSM_comb

figure(1);
%Compare to MSSM moment rate estimates
subplot(1,2,1)
plot([14.5 17.5],[14.5 17.5],'k--'); hold on
%plots for case where ss=1 and ss=3 (i.e char cases with different widths)
plot(log10(num_MSSM(:,13)),log10(MoRate_StochasticEventCatalog{1}(:,1)),'bx',...
log10(num_MSSM(:,14)),log10(MoRate_StochasticEventCatalog{3}(:,1)),'rx','MarkerSize',12,'LineWidth',1.5);

axis([14.5 17.5 14.5 17.5]); axis square; hold on;
xlabel('MSSM Moment Rate (log Nm/yr)'),ylabel('Stochastic Event Catalog Moment Rate (log Nm/yr)'); hold on;
%legend('y=x','length-limited width','layer-limited width','Location','southeast');
grid on; set(gca,'fontsize',13)

save('MSSM_adapted_comparison','MoRate_StochasticEventCatalog','num_MSSM');

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

%load data from MSSM fault-based syncat

all_mag_range_GR = ((Mmin-dm/2):0.05:8.5)'; bg_AnnualRate = zeros(length(all_mag_range_GR),1);

MSSM_adaptedAnnualRate=zeros(length(all_mag_range_GR),length(syncat_opt)); 


for ll = 1:length(all_mag_range_GR)
    bg_AnnualRate (ll,1) = length(find(bg_catalog(:,4) >= all_mag_range_GR(ll)))/(NumSimu*t_limit); %Rate in Malawi from bg catalog
    for ss=1:length(syncat_opt)
    MSSM_adaptedAnnualRate(ll,ss) = length(find(StochasticEventCatalog{ss}(:,4) >= all_mag_range_GR(ll)))/(NumSimu*t_limit);%combined annual rate of events
    MSSM_adaptedAnnualRate_comb(ll,ss) = MSSM_adaptedAnnualRate(ll,ss)+bg_AnnualRate(ll);
    end
    
    fault_syncatRate (ll,1) = length(find(EQCAT(:,5) >= all_mag_range_GR(ll)))/(NumSimu*t_limit);%Rate from MSSM direct catalog

end

%% Make MFD Plot comparison

figure(2);

%Comparing all different fault-based catalogs
semilogy(all_mag_range_GR,MSSM_adaptedAnnualRate(:,1),'k-','LineWidth',1.2); hold on; %MFD from event catalog
semilogy(all_mag_range_GR,MSSM_adaptedAnnualRate(:,2),'k--x','LineWidth',1.2); hold on; %MFD from event catalog
semilogy(all_mag_range_GR,MSSM_adaptedAnnualRate(:,3),'-','Color',[0.6 0.6 0.6],'LineWidth',1.2); hold on; %MFD from event catalog
semilogy(all_mag_range_GR,MSSM_adaptedAnnualRate(:,4),'--x','Color',[0.6 0.6 0.6],'LineWidth',1.2); hold on; %MFD from event catalog
semilogy(all_mag_range_GR,fault_syncatRate,'b-','LineWidth',1.2); hold on; %MFD from event catalog
semilogy(all_mag_range_GR,bg_AnnualRate,'r-','LineWidth',1.2); hold on; %MFD from event catalog
set(gca,'fontsize',13);  axis([5.5 8.5 10^-5 10^-1]); legend('length-limited w, char','length-limited w, G-R',...
    'layer-limited w, char','layer-limited w, GR','MSSM Direct','background','Location','southwest');hold on;
xlabel('Magnitude'); ylabel('Annual frequency'); grid on; axis square;

%% Moment rate comparisons between MSSM, YC-85 curves, and catalogs

MoRate_Comp=zeros(length(StochasticEventCatalog),3);
source_mo=zeros(length(YC85SeismicMoment_char),3);

%Get MoRate from catalogs in MSSM
MoRate_Comp(1:2,1)=ones(2,1).*sum(MSSM_adapt.Var13);
MoRate_Comp(3:4,1)=ones(2,1).*sum(MSSM_adapt.Var14);

%Get Morate from MSSM curves
for ii=1:length(YC85SeismicMoment_char)
    source_mo(ii,1:2)=sum(YC85SeismicMoment_char{ii}.*RecurrenceVar{ii}(:,5));
    source_mo(ii,3:4)=sum(YC85SeismicMoment_exp{ii}.*RecurrenceVar{ii}(:,5));
end

MoRate_Comp(1,2)=sum(source_mo(:,1));%char, length limited
MoRate_Comp(2,2)=sum(source_mo(:,3));%G-R, length limited
MoRate_Comp(3,2)=sum(source_mo(:,2));%char, layer limited
MoRate_Comp(4,2)=sum(source_mo(:,4));%char, layer limited

for ii=1:length(StochasticEventCatalog)
    
    %Get MoRate from catalogs
    tmp_Mo= 10.^(1.5*StochasticEventCatalog{ii}(:,4)+9.05);
    MoRate_Comp(ii,3)= sum(tmp_Mo)/(NumSimu*t_limit);

end


%% Freq mag plot for particular fault or multifault

%Plots both theoritical curves and actual curves for multiple faults

ff=[314 406]; %select fault or multi_fault based on ID
%Note these plots do not work for sources that have rupture weightings <1!!!!
ymax=[10^-3,10^-1];%these parameters may need to change depending on fault
ymin=[10^-6,10^-5];%these parameters may need to change depending on fault
plt_label=vertcat(["(a)","(b)"],["(c)","(d)"]);

fntsize=9;
legendfntsize=7;

figure(2000);
tiledlayout(2,2,'tilespacing','compact');

for gg=1:length(ff)
    
    nexttile  
      
    indx=find(num_MSSM(:,1)==ff(gg));

    flt_ruptures={}; sec_ruptures={};
    
    if ff(gg)>600
        %finf multifault ruptures in EQCAT
        m_flt_indx = find(ff(gg)==num_multi_fault(:,1));
        m_flt_ruptures=find(EQCAT(:,7)==ff(gg));
        
        %find participating fault of multifault system
        flt_id=num_fault(find(ff(gg)==num_fault(:,23)),1);
    
        %find fault ruptures in MSSM
        for i=1:length(flt_id)
            flt_indx(i)=find(num_fault(:,1)==flt_id(i));
            if num_fault(flt_indx(i),13)>0%if fault is segmented
                %source weightings
                weighting_f = weight_SWM(1,1);
                weighting_s = weight_SWM(1,2);
                weighting_m = weight_SWM(1,3);
            else
                weighting_f = weight_SWM(1,1)/(weight_SWM(1,1)+weight_SWM(1,3));
                weighting_m = weight_SWM(1,3)/(weight_SWM(1,1)+weight_SWM(1,3));
            end
        end
    
        %find sec in MSSM
        sec_id=num_sec1(find(ff(gg)==num_sec1(:,18)),1);  
        if isempty(sec_id)==0;
            for i=1:length(sec_id)
            %Index num_sec1 where no linking sections
            sec_indx(i)=find(num_sec1(:,1)==sec_id(i));
            end 
        end
    
        %Plot actual and catalog multifault ruptures 
        m1=semilogy(num_multi_fault(m_flt_indx,10),(1/num_multi_fault(m_flt_indx,13)),'kv','LineWidth',1.1); hold on; %Plot M-F from MSSM for multi_fault
        
        mcat_rate=length(find(EQCAT(:,7)==ff(gg)))/(NumSimu*weighting_m);
        m2=semilogy(num_multi_fault(m_flt_indx,10),mcat_rate,'r^','LineWidth',1.1); hold on;
        
        for i=1:length(flt_id)
            f1=semilogy(num_fault(flt_indx(i),18),(1/num_fault(flt_indx(i),21)),'kv','LineWidth',1.1); hold on; %Plot M-F from MSSM for fault
            fcat_rate=length(find(EQCAT(:,7)==flt_id(i)))/(NumSimu*weighting_f);
            f2=semilogy(num_fault(flt_indx(i),18),fcat_rate,'r^','LineWidth',1.1); hold on; %Plot M-F from MSSM for fault
        end
        
         if isempty(sec_id)==0 %i.e. fault is segmented
        
            for i=1:length(sec_id)
                s1(i)=semilogy(num_sec1(sec_indx(i),21),(1/num_sec1(sec_indx(i),24)),'kv','LineWidth',1.1); hold on; %Plot M-F from MSSM for sec
                scat_rate=length(find(EQCAT(:,7)==sec_id(i)))/(NumSimu*weighting_s);
                s2(i)=semilogy(num_sec1(sec_indx(i),21),scat_rate,'r^','LineWidth',1.1); hold on; %Plot M-F from MSSM for fault
            end
        end
    
    else %if not a multifault
    
        %find if participating sections in fault system
        flt_indx=find(ff(gg)==num_fault(:,1));
        if num_fault(flt_indx,13)>0
            sec_id=num_sec1(find(ff(gg)==num_sec1(:,15)),1); 
            for i=1:length(sec_id)
                %Index num_sec1 where no linking sections
                sec_indx(i)=find(num_sec1(:,1)==sec_id(i));
            end
            weighting_s = weight_SWM(1,2)/(weight_SWM(1,1)+weight_SWM(1,2));%weighting between fault and sec ruptures
            weighting_f = weight_SWM(1,1)/(weight_SWM(1,1)+weight_SWM(1,2));
        else
            weighting_f=1;
        end
    
     f1=semilogy(num_fault(flt_indx,18),(1/num_fault(flt_indx,21)),'kv','LineWidth',1.1); hold on; %Plot M-F from MSSM for fault
     fcat_rate=length(find(EQCAT(:,7)==ff(gg)))/(NumSimu*weighting_f);
     f2=semilogy(num_fault(flt_indx,18),fcat_rate,'r^','LineWidth',1.1); hold on; %Plot M-F from MSSM for fault
     
     if isempty(sec_id)==0 %i.e. fault is segmented
        for i=1:length(sec_id)
            s1(i)=semilogy(num_sec1(sec_indx(i),21),(1/num_sec1(sec_indx(i),24)),'kv','LineWidth',1.1); hold on; %Plot M-F from MSSM for sec
            scat_rate=length(find(EQCAT(:,7)==sec_id(i)))/(NumSimu*weighting_s);
            s2(i)=semilogy(num_sec1(sec_indx(i),21),scat_rate,'r^','LineWidth',1.1); hold on; %Plot M-F from MSSM for fault
        end

    end 

    end
    
    sourceRate={}; source_events={};

    %Index events from ff in StochasticEventCatalog
    for ss=1:length(StochasticEventCatalog)
        source_events{ss}=find(StochasticEventCatalog{ss}(:,6) == ff(gg)); 
        sourceRate{ss}=zeros(length(mag_range_GR{indx}),1);
         %Index events based on recurrence model and fault width
            for ii = 1:length(mag_range_GR{indx})
                sourceRate{ss}(ii,1) = length(find(StochasticEventCatalog{ss}(source_events{ss},4) >= mag_range_GR{indx}(ii)))/(NumSimu*t_limit);%annual rate of events where width=1
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

    %Plot for intermediate estimates
    %plot width =1
    a1=semilogy(mag_range_GR{indx},fm_char_var1{indx}(:,((length(RecurrenceVar{indx})/2)+0.5)*3),'b-','LineWidth',1.1); hold on;
    semilogy(mag_range_GR{indx},fm_exp_var1{indx}(:,((length(RecurrenceVar{indx})/2)+0.5)*3),'b-' ,'LineWidth',1.1); hold on; 
    %plot width =2
    semilogy(mag_range_GR{indx},fm_char_var2{indx}(:,((length(RecurrenceVar{indx})/2)+0.5)*3),'b-','LineWidth',1.1); hold on
    semilogy(mag_range_GR{indx},fm_exp_var2{indx}(:,((length(RecurrenceVar{indx})/2)+0.5)*3),'b-','LineWidth',1.1); hold on;
    
    %MFD from event catalog
    
    for ss=1:length(StochasticEventCatalog)
        a2=semilogy(mag_range_GR{indx},sourceRate{ss},'LineWidth',1.1,'color',[0.5 0.5 0.5]);  
    end
      
    legend([f1,f2,a1,a2],['Direct MSSM'],['Direct MSSM' newline 'Simulated Events'],...
        ['Adapted MSSM' newline 'Recurrence Models'],['Adapted MSSM' newline 'Simulated Events'],...
      'FontSize',legendfntsize,'Location','southwest'); hold on;
    set(gca,'fontsize',fntsize); axis([6 num_MSSM(indx,12)+0.2 ymin(gg) ymax(gg)]); hold on;
    xlabel('Magnitude'); ylabel('Annual frequency of exceedance'); grid on; axis square;

    t1 = title(plt_label(gg,1),'FontWeight','normal');
    set(t1, 'horizontalAlignment', 'left','units','normalized');
    h1 = get(t1, 'position');
    set(t1, 'position', [-0.2 h1(2) h1(3)]); set(gca,'FontSize',fntsize)

    nexttile
    %Plot for weighted average from recurrence models
    % Note weighted average cuts out at mmax-mshift (i,e lower bound of Mmax)
    %plot width =1
    a1=semilogy(mag_range_GR{indx},weighted_GR_char1,'b-','LineWidth',1.1); hold on
    semilogy(mag_range_GR{indx},weighted_GR_exp1,'b-' ,'LineWidth',1.1); hold on; 
    %plot width =2
    semilogy(mag_range_GR{indx},weighted_GR_char2,'b-','LineWidth',1.1); hold on
    semilogy(mag_range_GR{indx},weighted_GR_exp2,'b-','LineWidth',1.1); hold on;
    %MFD from event catalog
    for ss=1:length(StochasticEventCatalog)
        a2=semilogy(mag_range_GR{indx},sourceRate{ss},'LineWidth',1.1,'color',[0.5 0.5 0.5]);  
    
    end      
    legend([a1, a2],['Adapted MSSM' newline 'Recurrence Models'],['Adapted MSSM' newline 'Simulated Events'],...
       'FontSize',legendfntsize,'Location','southwest'); hold on;
    set(gca,'fontsize',fntsize); axis([6 num_MSSM(indx,12)+0.2 ymin(gg) ymax(gg)]); hold on;
    xlabel('Magnitude'); ylabel('Annual frequency of exceedance'); grid on; axis square;


    t1 = title(plt_label(gg,2),'FontWeight','normal');
        set(t1, 'horizontalAlignment', 'left','units','normalized');
    h1 = get(t1, 'position');
    set(t1, 'position', [-0.2 h1(2) h1(3)]); set(gca,'FontSize',fntsize)

end

set(gcf,'Position',  [183 263 590 534]);

%% Plot adapted_MSSM results with MSSM itself and direct MSSM eqcat catalog 
%Only considers fault width = 1 (for concistency with MSSM database)

gg=1; %only plots for one fault as selected from variable ff

figure(5);

indx=find(num_MSSM(:,1)==ff(gg));

a1=semilogy(mag_range_GR{indx},fm_char_var1{indx}(:,((length(RecurrenceVar{indx})/2)+0.5)*3),'r--',mag_range_GR{indx},fm_exp_var1{indx}(:,((length(RecurrenceVar{indx})/2)+0.5)*3),'k--' ,'LineWidth',1.5); hold on; 
%plot width =2
a2=semilogy(mag_range_GR{indx},fm_char_var2{indx}(:,((length(RecurrenceVar{indx})/2)+0.5)*3),'r-',mag_range_GR{indx},fm_exp_var2{indx}(:,((length(RecurrenceVar{indx})/2)+0.5)*3),'k-','LineWidth',1.5); hold on;

if ff(gg)>600 %For multi_fault, plot results from MSSM and eqcat
     
    m_flt_indx = find(ff(gg)==num_multi_fault(:,1));
    
    m1=semilogy(num_multi_fault(m_flt_indx,10),(1/num_multi_fault(m_flt_indx,13)),'k*','Markersize',8,'LineWidth',1.3); hold on; %Plot M-F from MSSM for multi_fault
    
    flt_id=num_fault(find(ff(gg)==num_fault(:,23)),1);
    
    for i=1:length(flt_id)
        flt_indx(i)=find(num_fault(:,1)==flt_id(i));
        f1(i)=semilogy(num_fault(flt_indx(i),18),(1/num_fault(flt_indx(i),21)),'b*','LineWidth',1.3,'Markersize',8); hold on; %Plot M-F from MSSM for fault   
      if i>1 %avoid duplication for legend indexing
         set(get(get(f1(i),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
      end
    end
   
   sec_id=num_sec1(find(ff(gg)==num_sec1(:,18)),1);  
   if isempty(sec_id)==0;
           for i=1:length(sec_id)
            %Index num_sec1 where no linking sections
            sec_indx(i)=find(num_sec1(:,1)==sec_id(i));
           end 

        
    for i=1:length(sec_id)
        s1(i)=semilogy(num_sec1(sec_indx(i),21),(1/num_sec1(sec_indx(i),24)),'*','LineWidth',1.3,'Markersize',8,'MarkerEdgeColor',[0.4 0.4 0.4]); hold on; %Plot M-F from MSSM for sec
      if i>1  %avoid duplication for legend indexing
         set(get(get(s1(i),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
      end
        
    end
        legend([s1(1) f1(1) m1 a1(1) a1(2) a2(1) a2(2)],{['MSSM MultiFault'], ['MSSM Fault'],['MSSM Section'],['Adapted char events:' newline 'length-limited width'],['Adapted G-R events:' newline 'length-limited width'],...
            ['Adapted char events:' newline 'layer-limited width'],['adapted G-R events:' newline 'layer-limited width']},...
            'Location','northeast','fontsize',12); hold on;
    else
         legend([m1 f1 a1(1) a1(2) a2(1) a2(2)],{['MSSM MultiFault'], ['MSSM Fault'],['Adapted char events:' newline 'length-limited width'],['Adapted G-R events:' newline 'length-limited width'],...
             ['Adapted char events:' newline 'later-limited width'],['Adapted G-R events:' newline 'layer-limited width']},...
            'Location','northeast','fontsize',12); hold on;
    end
    
else %fault is not part of multi_fault system
    
    flt_indx=find(ff(gg)==num_fault(:,1));
    f1=semilogy(num_fault(flt_indx,18),(1/num_fault(flt_indx,21)),'b*','LineWidth',1.3,'Markersize',8); hold on; %Plot M-F from MSSM for fault   
      
    sec_id=num_sec1(find(ff(gg)==num_sec1(:,15)),1);  
    
    if isempty(sec_id)==0
        for i=1:length(sec_id)
            s1=semilogy(num_sec1(sec_indx(i),21),(1/num_sec1(sec_indx(i),24)),'*','LineWidth',1.3,'Markersize',8,'MarkerEdgeColor',[0.4 0.4 0.4]); hold on; %Plot M-F from MSSM for sec
        
            legend([s1 f1 a1(1) a1(2) a2(1) a2(2)],{['MSSM Section'],['MSSM Fault'],['Adapted char events:' newline 'length-limited width'],['Adapted G-R events:' newline 'length-limited width'],...
            ['Adapted char events:' newline 'layer-limited width'],['Adapted G-R events:' newline 'layer-limited width']},...
             'Location','northeast','fontsize',11); hold on;
        end
    else
        legend([f1 a1(1) a1(2) a2(1) a2(2)],{['MSSM Fault'],['Adapted char events:' newline 'length-limited width'],['Adapted G-R events:' newline 'length-limited width'],...
            ['Adapted char events:' newline 'layer-limited width'],['Adapted G-R events:' newline 'layer-limited width']},...
            'Location','northeast','fontsize',11); hold on;
    end 

end

axis([5 num_MSSM(indx,12)+0.2 10^-6 10^-2]); xlabel('Magnitude'); ylabel('Annual occurrence rate'); axis square; grid on; set(gca,'fontsize',15)

%% Plot all theortical model curves for the MSSM source and MO rate comparison

figure(6);
gg=1;
indx=find(num_MSSM(:,1)==ff(gg));
fntsize=10;

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

tiledlayout(1,3);
nexttile

p1=semilogy(mag_range_GR{indx},fm_char_var1{indx}(:,3:3:end),'Color',[0.5 0.5 0.5],'LineWidth',1); hold on
p2=semilogy(mag_range_GR{indx},weighted_GR_char1,'k-','LineWidth',1.5); hold on; 

xlabel('Magnitude'); ylabel('Annual frequency of exceedance');
legend([p1(1) p2],'Individual','Weigted Average','Location','Southwest');
axis([6 num_MSSM(indx,12)+0.2 10^-6 10^-3]); axis square; grid on
t1 = title('(a)','FontWeight','normal');
set(t1, 'horizontalAlignment', 'left','units','normalized');
h1 = get(t1, 'position');
set(t1, 'position', [-0.2 h1(2) h1(3)]); set(gca,'FontSize',fntsize)

nexttile

p1=semilogy(mag_range_GR{indx},fm_exp_var1{indx}(:,3:3:end),'Color',[0.5 0.5 0.5],'LineWidth',1); hold on
p2=semilogy(mag_range_GR{indx},weighted_GR_exp1,'k-','LineWidth',1.5); hold on; 

xlabel('Magnitude'); ylabel('Annual frequency of exceedance');
legend([p1(1) p2],'Individual','Weigted Average','Location','Southwest');
axis([6 num_MSSM(indx,12)+0.2 10^-6 10^-3]); axis square; grid on
t1 = title('(b)','FontWeight','normal');
set(t1, 'horizontalAlignment', 'left','units','normalized');
h1 = get(t1, 'position');
set(t1, 'position', [-0.2 h1(2) h1(3)]); set(gca,'FontSize',fntsize)

nexttile %plot theoritical seismic moment rates against each other

tmp_indx=find(MSSM_adapt.Var1==ff(gg));
target_mo(1)=MSSM_adapt.Var13(tmp_indx); %index target Mo for length limited
target_mo(2)=MSSM_adapt.Var14(tmp_indx);%index target Mo for layer limited

p1=plot([1e15 1e16],[1e15 1e16],'k--','linewidth',1.05); hold on
p2=plot(target_mo(1),target_mo(1),'k*','linewidth',1.5,'markersize',8);
plot(target_mo(2),target_mo(2),'k*','linewidth',1.5,'markersize',8);
p4=plot(YC85SeismicMoment_char{indx}(:,1),YC85SeismicMoment_exp{indx}(:,1),'bv','linewidth',1.1);
p5=plot(YC85SeismicMoment_char{indx}(:,2),YC85SeismicMoment_exp{indx}(:,2),'r^','linewidth',1.1); axis square
xlabel({'Characteristic model', 'moment rate (Nm/yr)'}); ylabel('G-R model moment rate (Nm/yr)'); axis([2e15 5e15 2e15 5e15]); 
grid on
legend([p1,p2,p4,p5],{'y=x','Target Moment Rate','Length limited width','Layer limited width'},'location','northwest');

t1 = title('(c)','FontWeight','normal');
set(t1, 'horizontalAlignment', 'left','units','normalized');
h1 = get(t1, 'position');
set(t1, 'position', [-0.2 h1(2) h1(3)]); set(gca,'FontSize',fntsize)

set(gcf,'Position',[440 467 911 330]);

MoRate_Comp_flt=horzcat(RecurrenceVar{indx}(:,3:5),YC85SeismicMoment_char{indx}(:,1)/1e15,YC85SeismicMoment_exp{indx}(:,1)/1e15,(YC85SeismicMoment_exp{indx}(:,1)./YC85SeismicMoment_char{indx}(:,1)));