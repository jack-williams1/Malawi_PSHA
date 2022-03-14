%% MSSD Direct Rupture Weighting

clear
close all

%add in other paths
mydir  = pwd; idcs   = strfind(mydir,'/');
addpath(mydir(1:idcs(end)-1));
addpath([mydir(1:idcs(end)-1) '/misc_functions']);

% Finite fault model
load MSSD_sources

load('syncat_PSHA_MSSD_input','para_Mw','NumSimu','t_limit');


%Remove linking sections
%i.e those sections who don't rupture on their own and have a moment rate set as '0'
%save in new variables num_fault1 and num_sec1

tmp_bb=find(num_sec(:,26)); num_sec1=num_sec(tmp_bb,:); clear tmp_bb
tmp_bb=find(num_fault(:,25)); num_fault1=num_fault(tmp_bb,:);

%% Generate exhaustive list of all possible rupture weightings 

weight_increment = 0.1;
clear weight_SFM_comb
% Weight related to the rupture type: (1) fault, (2) section, (3) multifault 

num_SFM_comb = 0;

    for ii1 = 1:1/weight_increment+1
        for ii2 = 1:1/weight_increment+1
            for ii3 = 1:1/weight_increment+1
                tmp = [(ii1-1)*weight_increment (ii2-1)*weight_increment (ii3-1)*weight_increment];
                %Finds all times when 3 weighting components that add up to
                %1 and weighting does not equal -
                if sum(tmp) == 1 && (ii1-1)*weight_increment>0 && (ii2-1)*weight_increment>0 && (ii3-1)*weight_increment>0  
                    num_SFM_comb = num_SFM_comb + 1;
                    weight_SFM_comb(num_SFM_comb,1:3) = tmp;  %Saves combinations when weighting components add up to 1
                end
            end
       
        end
    end
        
weight_SFM_comb=vertcat(weight_SFM_comb,[weight_SFM_comb(:,1) weight_SFM_comb(:,3), weight_SFM_comb(:,2)]); 
weight_SFM_comb=unique(weight_SFM_comb,'rows');

%% Generate 5 recurrence models

%recurrence_para1: mag and RI for MSSD faults
%recurrence_para2: mag and RI for MSSD sections and MSSD faults if faults unsegmented
%recurrence_para3: mag and RI for MSSD sections and MSSD multifaults and faults if unsegmented
%recurrence_para4: mag and RI for MSSD multifauts and MSSD faults if faults isolated
%recurrence_para5: mag and RI for MSSD multifaults and MSSD sections if faults segmented and isolated

%recurrence_para1: index fault int mag, RI, and ID
recurrence_para1=[num_fault1(:,18), num_fault1(:,21)./num_fault1(:,3) num_fault1(:,1)];

%RECURRENCE PARA FOR SECTION RUPTURES

f_s_match={};
m_s_match={};
f_m_match={};

%find if fault is unsegmented (i.e is also represented in sec catalog)
    for ff=1:length(num_fault1)
        f_s_match{ff}=find(num_fault1(ff,1)==num_sec1(:,15)); 
   
    if isempty(f_s_match{ff})==1
        recurrence_para_tmp2(ff,:)=[num_fault1(ff,18), num_fault1(ff,21)./num_fault1(ff,3), num_fault1(ff,1)];
    end
 
    end
    
recurrence_para_tmp2( ~any(recurrence_para_tmp2,2), : ) = [];

recurrence_para_tmp3=zeros(1,3); tmpindx=1;    
    
%find if multi_fault contains discrete fault sections
 for mm=1:length(num_multi_fault)
        m_s_match{mm}=find(num_multi_fault(mm,1)==num_sec1(:,18)); 
   
    if isempty(m_s_match{mm})==1
        %Include multi_fault ruptures that don't contain discrete fault
        %sections
        recurrence_para_tmp3(tmpindx,:)=[num_multi_fault(mm,10),num_multi_fault(mm,13)./num_multi_fault(mm,3),num_multi_fault(mm,1)];
        tmpindx=tmpindx+1;
    end
 end
 
%Include ruptures that are neither segmented faults or multi faults AND
%include faults of multi fault systems, in which other faults are segmented 
 
recurrence_para_tmp4=zeros(1,3); tmpindx=1;
  for ff=1:length(recurrence_para_tmp2)
        tmp=find(recurrence_para_tmp2(ff,3)==num_fault(:,1)); 
        
        if num_fault(tmp,23)==0
            recurrence_para_tmp4(tmpindx,:)=recurrence_para_tmp2(ff,:);
            tmpindx=tmpindx+1;
        end
  end  
  
recurrence_para2=vertcat(recurrence_para_tmp2(:,1:3),[num_sec1(:,21),num_sec1(:,24)./num_sec1(:,5),num_sec1(:,1)]);
recurrence_para3=vertcat(recurrence_para_tmp3(:,1:3),recurrence_para_tmp4(:,1:3),[num_sec1(:,21),num_sec1(:,24)./num_sec1(:,5),num_sec1(:,1)]);

clear recurrence_para_tmp2 recurrence_para_tmp3 recurrence_para_tmp4
clear m_s_match m_f_s_match

%RECURRENCE PARA FOR MULTI-FAULT RUPTURES

f_m_match={}; m_s_match={};

recurrence_para_tmp4=zeros(1,3); 
tmpindx=1;

%find if fault part of multi-fault system
    for ff=1:length(num_fault1)
        f_m_match{ff}=find(num_fault1(ff,23)==num_multi_fault(:,1)); 
%select faults not part of multi-fault system   
    if isempty(f_m_match{ff})==1
        recurrence_para_tmp4(tmpindx,:)=[num_fault(ff,18), num_fault(ff,21)./num_fault(ff,3), num_fault(ff,1)];
        tmpindx=tmpindx+1;
    end
 
    end
recurrence_para_tmp5=zeros(1,3); 
recurrence_para_tmp6=zeros(1,3); 

tmpindx1=1; tmpindx2=1; tmp_idx={};

%find if non-multifault fault is segmented
for ss=1:length(recurrence_para_tmp4)
        m_s_match{ss}=find(recurrence_para_tmp4(ss,3)==num_sec1(:,15)); 
   
    if isempty(m_s_match{ss})==1
        recurrence_para_tmp5(tmpindx1,:)=recurrence_para_tmp4(ss,:);%keep with fault id false
        tmpindx1=tmpindx1+1;
    else
        tmp_idx{ss}=find(recurrence_para_tmp4(ss,3)==num_sec1(:,15));
            for kk=1:height(tmp_idx{ss})
                 recurrence_para_tmp6(tmpindx2,:)=[num_sec1(tmp_idx{ss}(kk),21),num_sec1(tmp_idx{ss}(kk),24)./num_sec1(tmp_idx{ss}(kk),5),num_sec1(tmp_idx{ss}(kk),1)];%take section id if true
                 tmpindx2=tmpindx2+1;
            end  
    end
end  

recurrence_para4=vertcat(recurrence_para_tmp4(:,1:3),[num_multi_fault(:,10),num_multi_fault(:,13)./num_multi_fault(:,3),num_multi_fault(:,1)]);
recurrence_para5=vertcat(recurrence_para_tmp5(:,1:3),recurrence_para_tmp6,[num_multi_fault(:,10),num_multi_fault(:,13)./num_multi_fault(:,3),num_multi_fault(:,1)]);

clear recurrence_para_tmp5 recurrence_para_tmp6

%% Run syncat simulations through each rupture weighting combination

for ijk=1:length(weight_SFM_comb)

clear weight_SWM

simu_count = 1;

global_count = 0; 

clock_old = cputime;

Num_event = zeros(NumSimu,1);

Count_WS  = zeros(3,1);    
    
%loop through each possible combination of rupture types
weight_SWM(1,1:3) = weight_SFM_comb(ijk,:);   
weight_SWM(2,1:3) = 1:3;

for ii = 1:NumSimu
    
    if ii == 400000*simu_count
        clock_new = cputime;
        disp(['Simulation cycle is: ',num2str(ii),' of combination ',num2str(ijk),' & Global count = ',num2str(global_count),' & Required time is: ',num2str(clock_new-clock_old)]);
        clock_old = clock_new;
        simu_count = simu_count + 1;
    end
    
    % Select a rupture regime for each simulation cycle
        ran = rand(1,1);
        for jj = 1:3
            if  ran > sum(weight_SWM(1,1:jj-1)) && ran <= sum(weight_SWM(1,1:jj)) %select weighting based on value of ran
                  %and comparing it to value of rupture weightings
                  ws_pick      = weight_SWM(2,jj);
                  Count_WS(jj) = Count_WS(jj) + 1;
                  break
            end
        end
    
    % Magnitude and recurrence interval corresponding to the selected
    % rupture type (section, fault, or multi-fault)
        if ws_pick == 1
            recurrence_para = recurrence_para1; %Chose the magntiude and RI from whole fault ruptures list
        
        elseif ws_pick == 2 %chose magnitude and RI from section rupture list
                
        %for section ruptures, chose weighting between faults and
        %multifaults consistent with weight_SWM 
                tmp=sum(weight_SWM(1,1)+weight_SWM(1,3));
                ran_s = tmp*rand(1,1);  
                    if  ran_s < weight_SWM(1,1)
                        recurrence_para = recurrence_para2; %Chose parameters for unsegmented faults that are not multi faults
                    else
                        recurrence_para = recurrence_para3; %Chose parameters for unsegmented faults that are multi faults
                clear tmp
                    end
                    
        elseif ws_pick ==3 %chose magnitude and RI from multifault rupture list
            
                tmp=sum(weight_SWM(1,1)+weight_SWM(1,2));
                ran_m = tmp*rand(1,1);  
                    if  ran_m < weight_SWM(1,1)
                         recurrence_para = recurrence_para4; %Chose parameters for non mulitfaults that are faults 
                    else
                        recurrence_para = recurrence_para5; %Chose parameters for non mulitfaults that are also segmented
                    end
                clear tmp
        end

    count = 0;
    t_sum = 0;
    
    count_pre = global_count;

    occurrence_time = zeros(1,length(recurrence_para));%Creates 1 row of zeros with length equal to number of faults/sections
    
    while t_sum <= t_limit %During ii loop, create new set of t_sum values, only continue running if less than t_limit
        
        % Earthquake occurrence - Poisson process
        if count == 0
            occurrence_time = expinv(rand(1,length(recurrence_para)),recurrence_para(:,2)');
            %first creates a matrix of uniformly distributed random numbers between 0-1, where the number of columns is
            %equal to the number of faults or sections. i.e. represents a standard exponential distribution
            %For each random number, expinv finds the value x from an inverse cumulative distribution function where x=F^-1 (p/mu) and 
            % p is the random number, and mu is the earthquake recurrence interval, and represents the average value in the cdf
            %This is repeated for each fault/section, each of which is a different row in the matrix
    
        else
            occurrence_time(1,num) = occurrence_time(1,num) + expinv(rand(1,1),recurrence_para(num,2));
            % if count = 1, occurrence time in the first row equals 0 (as the first if statement will not be evaluated),
            %plus a random value drawn from the faults/section recurrence interval cumulative distribution
        end
        
        count = count + 1;
        
        % Find the smallest values among the occurrence times that were
        % calculated in the previous step for each column in the first row
        [t,num] = min(occurrence_time(1,:));
        t_sum = t; %if t_sum>t_limit, then while loop ends here
                    
        if t_sum <= t_limit %If one of those occurrence times is less than 1, record earthquake 
            
            % Earthquake magnitude - truncated normal distribution
            m = recurrence_para(num,1) + para_Mw(1)*randn(1,1); %select mag from fault or section and consider distribution
            m = min(max(m,recurrence_para(num,1)+para_Mw(2)),recurrence_para(num,1)+para_Mw(3));
                   
            rupture_id  = recurrence_para(num,3);  %Identify rupture id 
            global_count = global_count + 1;
                                                                      
            % Synthetic earthquake catalog
            % 1 ) Overall cycle
            % 2 ) Simulation cycle
            % 3 ) Earthquake magnitude
            % 4 ) Source ID 
            % 5 ) Section, whole fault, or multi fault rupture (1 to 3)
            EQCAT(global_count,1:5) = [global_count ii m rupture_id ws_pick];
            
        end %end if t_sum <= t_limit
    end %end while loop
    
    Num_event(ii,1) = global_count - count_pre; %Num_event for each year in catalog
    
end %end syncat simulation for each individual rupture weighting combination

EQCAT_comb{ijk}=EQCAT;


clear EQCAT

end

%% Save for later analysis

save('EQCAT_comb','EQCAT_comb','num_SFM_comb','weight_SFM_comb','NumSimu','t_limit');


%% Plot MFD of all curves

%load EQCAT_comb; uncomment if run previously
mydir  = pwd; idcs   = strfind(mydir,'/');
addpath(mydir(1:idcs(end)-1)); addpath([mydir(1:idcs(end)-1) '/syncat_bg']);
load syncat_bg

rem_event=find(bg_catalog(:,1)<0);
bg_catalog(rem_event,:)=[];

%% Regional seismicity (Zone 9 - Rukwa-Malawi Rift; Poggi et al. 2017) 
a_background = 4.7; % GEM - Poggi et al. (zone 9 scaled to Malawi)
b_background = 1.02;
mag_check = 4.5:0.05:8.5;

% Range over which MSSD is considered complete
l_cutoff =6; u_cutoff =7.6; bg_AnnualRate=zeros(length(mag_check),1);

for ijk = 1:length(weight_SFM_comb)
    for ii = 1:length(mag_check)
      Rate{ijk}(ii,1) = length(find(EQCAT_comb{ijk}(:,3)>=mag_check(ii)))/(NumSimu*t_limit); %Bin events cumulatively by magnitude at intervals of 0.05 
      if ijk==1
        bg_AnnualRate (ii,1) = length(find(bg_catalog(:,4) >= mag_check(ii)))/(NumSimu*t_limit);
      end
    end
end

figure(1);

%Plot all MFD curves
for ijk = 1:length(weight_SFM_comb)
    semilogy(mag_check,Rate{ijk}(:,1),'LineWidth',0.75); hold on;  
 end

%Highlight MFD curves for extreme scenarios
for i=1:length(weight_SFM_comb)
    
    if any(weight_SFM_comb(i,:,:)==0.8)
        p2=semilogy(mag_check,Rate{i}(:,1),'k-','LineWidth',1.8); hold on;  
    end
    
end

%Plot GEM background seismicity rate
p3=semilogy(mag_check,bg_AnnualRate,'r--','LineWidth',1.4);
line([l_cutoff,l_cutoff],[10^-5 10^0],'LineStyle','--','Color',[0.3 0.3 0.3],'LineWidth',1.2);hold on
line([u_cutoff,u_cutoff],[10^-5 10^0],'LineStyle','--','Color',[0.3 0.3 0.3],'LineWidth',1.2);
xlabel('Magnitude'); ylabel('Annual frequency of exceedance');
axis square; grid on; axis([5 9 10^-5 10^0]);
legend([p2 p3],{'MSSD Direct Extreme weightings','Areal Sources Simulated Catalog'});
set(gca,'FontSize',13.5);

%% Perform KS test on each curve with respect to theoritical curve

a_={};
for ijk = 1:length(weight_SFM_comb)
    
    %sample magnitudes only in range where catalog considered complete
    sample_mag{ijk}=EQCAT_comb{ijk}(find(EQCAT_comb{ijk}(:,3)>l_cutoff & EQCAT_comb{ijk}(:,3)<u_cutoff),3);
    
    ar=length(find(EQCAT_comb{ijk}(:,3)>l_cutoff))/NumSimu;%annual rate of events M>6
    a_{ijk} = log10(ar)+b_background*l_cutoff;%a value for record ijk
    mag_check_tmp = l_cutoff:0.01:u_cutoff;
    
   for mm=1:length(mag_check_tmp)
       %determine theoritical number of events in mag bin mm
       Rate_t{ijk}(mm)=(a_{ijk}-b_background*mag_check_tmp(mm));
       %determine simulated number of events in mag bin mm
       Rate_sim{ijk}(mm,1) = length(find(EQCAT_comb{ijk}(:,3)>=mag_check_tmp(mm)))/(NumSimu*t_limit); %
   end
   %ks test between simulated and theoritical events
   [h(ijk),p(ijk),ksdist(ijk)]=kstest2(Rate_sim{ijk},10.^Rate_t{ijk},'Alpha',0.1);
   
   clear ar

end

r_weight_table=table(weight_SFM_comb(:,1),weight_SFM_comb(:,2),weight_SFM_comb(:,3),round(ksdist',3),round(p',3),h','VariableNames',{'fault w','section w','multifault w','k-s dist','p','reject?'}) ;
writetable(r_weight_table,'r_weight_results.csv');

%% Plot min, intermediate, and max values of ksdist
[~,tmp]=min(abs(ksdist(:)-mean(ksdist)));
ks_indx=[find(ksdist==min(ksdist),1),tmp,find(ksdist==max(ksdist),1)];

figure(201);

for i=1:length(ks_indx)
    subplot(1,3,i)
    semilogy(mag_check_tmp,Rate_sim{ks_indx(i)},'k-','LineWidth',1); hold on
    semilogy(mag_check_tmp,10.^Rate_t{ks_indx(i)},'b-','LineWidth',1);  
    xlim([l_cutoff u_cutoff]);

    xlabel('Magnitude'); ylabel('Annual frequency of exceedance');
    axis square; grid on; axis([l_cutoff u_cutoff  10^-3 10^-1]);
    set(gca,'FontSize',13.5);
    
    legend('Simulated Trend','Theortical G-R');
    title({'Rupture Weight: ',num2str(weight_SFM_comb(ks_indx(i),:))});
    subtitle({'k-s distance =',num2str(ksdist(ks_indx(i)),3)});

end
set(gcf,'Position', [614 429 1240 454]);

%% Figure for manuscript that combines MFD and K-S tests

figure(301);
subplot(1,2,1)

ks_tmp1=max(find(ksdist==min(ksdist)));
ks_tmp2=24;%index for 0.4 0.3 0.3 weighting

%Plot all MFD curves
for ijk = 1:length(weight_SFM_comb)
 
        pp1=semilogy(mag_check,Rate{ijk}(:,1),'LineWidth',0.7,'Color',[0.6 0.6 0.6]); hold on;  
end


%{
Option to highlight MFD curves for extreme weighting combinations
for i=1:length(weight_SFM_comb)
    
    if any(weight_SFM_comb(i,:,:)==0.8)
        p2=semilogy(mag_check,Rate{i}(:,1),'b-','LineWidth',1); hold on;  
    end
    
end
%}

%Plot record with smallest k-s distance
pp2=semilogy(mag_check,Rate{ks_tmp1}(:,1),'k-','LineWidth',2); hold on;

l1=line([l_cutoff,l_cutoff],[10^-5 10^-1],'LineStyle','--','Color',[0.3 0.3 0.3],'LineWidth',1.2);hold on
line([u_cutoff,u_cutoff],[10^-5 10^-1],'LineStyle','--','Color',[0.3 0.3 0.3],'LineWidth',1.2);
xlabel('Magnitude'); ylabel('Annual frequency of exceedance');
legend([pp1 pp2 l1],'all weightings','selected weighting',['assessed magnitude' newline 'range'])
axis square; grid on; axis([5 8.5 10^-5 10^-1]);
set(gca,'FontSize',13.5);

t1 = title('(a)','fontweight','normal','FontSize',15);
set(t1, 'horizontalAlignment', 'left','units','normalized');
h1 = get(t1, 'position');
set(t1, 'position', [-0.2 h1(2) h1(3)]);

ks_indx=[ks_tmp1,ks_tmp2];
col_opt=['b','r'];

%Detailed plot with smallest K-S distance and base-case of equal weightings
subplot(1,2,2)

for i=1:2
    dummy(i)=semilogy(mag_check_tmp,nan(size(mag_check_tmp))); hold on %create dummy variable for legend
    p2(i)=semilogy(mag_check_tmp,Rate_sim{ks_indx(i)},'-','Color',col_opt(i),'LineWidth',1.2); hold on
    p3(i)=semilogy(mag_check_tmp,10.^Rate_t{ks_indx(i)},'--','Color',col_opt(i),'LineWidth',1.2);  
    xlim([l_cutoff u_cutoff])
   
end

xlabel('Magnitude'); ylabel('Annual frequency of exceedance');
axis square; grid on; axis([l_cutoff u_cutoff  10^-3 10^-1]);
set(gca,'FontSize',13.5);

[l1, objH]=legend({[newline 'junk1'],'simulated catalog','b-value =1.02',[newline 'junk2'],'simulated catalog','b-value =1.02'});

% Edit legend so contains subtitles with details of rupture weightings and ks-dist
set(findobj(objH, 'Tag', [newline 'junk1']), 'Vis', 'off');  set(findobj(objH, 'Tag', [newline 'junk2']), 'Vis', 'off');           % Make "junk" lines invisible
pos1 = get(objH(1), 'Pos');  pos2 = get(objH(4), 'Pos');                    
set(objH(1), 'Pos', [0.05 pos1(2:3)], 'String', ['Rupture Weight: ' num2str(weight_SFM_comb(ks_indx(1),2)) ':' num2str(weight_SFM_comb(ks_indx(1),1)),...
    ':' num2str(weight_SFM_comb(ks_indx(1),3)) newline 'K-S distance = ',num2str(ksdist(ks_indx(1)),2)]);
set(objH(4), 'Pos', [0.05 pos2(2:3)], 'String', ['Rupture Weight: ' num2str(weight_SFM_comb(ks_indx(2),2)) ':' num2str(weight_SFM_comb(ks_indx(2),1)),...
    ':' num2str(weight_SFM_comb(ks_indx(2),3)) newline 'K-S distance = ',num2str(ksdist(ks_indx(2)),2)]);
l1.Position= [0.745 0.6219 0.15 0.2856];

t1 = title('(b)','Fontweight','normal','FontSize',15);
set(t1, 'horizontalAlignment', 'left','units','normalized');
h1 = get(t1, 'position');
set(t1, 'position', [-0.2 h1(2) h1(3)]);


 set(gcf,'Position',[681 536 1100 443]);