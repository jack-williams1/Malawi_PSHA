%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%   Synthetic Catalog Generation for Fault Sources in Malawi   %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%For direct interpretation of the MSSM%
%ie directly uses estimates of recurrence intervals and magnitudes in the
%MSSM
%Considers single section, fault, and multi fault ruptures

%Saves catalog in variable 'eqcat_test'

clear
close all

%add in other paths
mydir  = pwd; idcs   = strfind(mydir,'/');
addpath(mydir(1:idcs(end)-1));
addpath([mydir(1:idcs(end)-1) '/misc_functions']);
addpath([mydir(1:idcs(end)-1) '/gis_files']);

load ('syncat_PSHA_MSSM_input','t_limit','NumSimu','syncat_name','weight_SWM','num_TD','weight_TD','para_Mw')

% Finite fault model
load MSSM_sources

%Remove linking sections
%i.e those sections who don't rupture on their own and have a moment rate set as '0'
%save in new variables num_fault and num_sec1
tmp_bb=find(num_sec(:,26)); num_sec1=num_sec(tmp_bb,:); clear tmp_bb

% Syncat inputs

%weight_SWM: Weight for (1) section rupture, (2) fault rupture, (3) multi-fault rupture
%para_Mw: Magnitude uncertainty: magnitude is normally distributed with (1) standard deviation, (2) lower truncation, and (3) upper truncation
%num_TD & weight_TD: % Depth uncertainty and weightings: multiple discrete values
%t_limit & NumSimu: Simulation length options

%% Generate 5 recurrence models

%recurrence_para1: int RI for MSSM faults
%recurrence_para2: int RI for MSSM sections and MSSM faults if faults unsegmented
%recurrence_para3: int RI for MSSM sections and MSSM multifaults and faults if unsegmented
%recurrence_para4: int RI for MSSM multifauts and MSSM faults if faults isolated
%recurrence_para5: int RI for MSSM multifaults and MSSM sections if faults segmented and isolated

%RECURRENCE PARA FOR FAULT SOURCES
%divide RI by weighting factor for ruptures on branching faults

recurrence_para1=[num_fault(:,18), num_fault(:,21)./num_fault(:,3) num_fault(:,1)];


%RECURRENCE PARA FOR SECTION SOURCES

f_s_match={};
m_s_match={};
f_m_match={};

tmpindx=1;

recurrence_para_tmp2=zeros(1,3);

%find if fault is unsegmented (i.e is also represented in sec catalog)
    for ff=1:length(num_fault)
        f_s_match{ff}=find(num_fault(ff,1)==num_sec1(:,15)); 
   
    if isempty(f_s_match{ff})==1
        recurrence_para_tmp2(tmpindx,:)=[num_fault(ff,18), num_fault(ff,21)./num_fault(ff,3), num_fault(ff,1)];
        tmpindx=tmpindx+1;
    end
 
    end

tmpindx=1;  

recurrence_para_tmp3=zeros(1,3); 
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
 
tmpindx=1;

recurrence_para_tmp4=zeros(1,3); 
%Include sources that are neither segmented faults or multi faults AND
%include faults of multi fault systems, in which other faults are segmented 
   
  for ff=1:length(recurrence_para_tmp2)
        tmp=find(recurrence_para_tmp2(ff,3)==num_fault(:,1)); 
        
        if num_fault(tmp,23)==0
            recurrence_para_tmp4(tmpindx,:)=recurrence_para_tmp2(ff,:);
            tmpindx=tmpindx+1;
        end
        clear tmp
  end

recurrence_para2=vertcat(recurrence_para_tmp2,[num_sec1(:,21),num_sec1(:,24)./num_sec1(:,5),num_sec1(:,1)]);

recurrence_para3=vertcat(recurrence_para_tmp3,recurrence_para_tmp4,...
    [num_sec1(:,21),num_sec1(:,24)./num_sec1(:,5),num_sec1(:,1)]);

clear recurrence_para_tmp3 recurrence_para_tmp4 recurrence_para_tmp5
clear m_s_match m_f_s_match

%RECURRENCE PARA FOR MULTI-FAULT SOURCES

f_m_match={}; m_s_match={};

recurrence_para_tmp4=zeros(1,3);
tmpindx=1;

%find if fault part of multi-fault system
 for ff=1:length(num_fault)
        f_m_match{ff}=find(num_fault(ff,23)==num_multi_fault(:,1)); 
%select faults not part of multi-fault system   
    if isempty(f_m_match{ff})==1
        recurrence_para_tmp4(tmpindx,:)=[num_fault(ff,18), num_fault(ff,21)./num_fault(ff,3), num_fault(ff,1)];
        tmpindx=tmpindx+1;
    end
 
 end

recurrence_para_tmp5=zeros(1,3); recurrence_para_tmp6=zeros(1,3);
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

recurrence_para4=vertcat(recurrence_para_tmp4,[num_multi_fault(:,10),num_multi_fault(:,13)./num_multi_fault(:,3),num_multi_fault(:,1)]);

recurrence_para5=vertcat(recurrence_para_tmp5,recurrence_para_tmp6,...
    [num_multi_fault(:,10),num_multi_fault(:,13)./num_multi_fault(:,3),num_multi_fault(:,1)]);
   
clear recurrence_para_tmp4 recurrence_para_tmp5 recurrence_para_tmp6
clear m_s_match m_f_s_match



%% Synthetic finite fault earthquake catalog

%[~,save_file] = xlsread(inputexcel,1,'C20:C20');

simu_count = 1;
global_count = 0; 
clock_old = cputime;
Num_event = zeros(NumSimu,1);
Count_WS  = zeros(3,1);

Count_EventW=0;
Count_EventS_1 =0; Count_EventS_2 =0; 
Count_EventM_1 =0; Count_EventM_2 =0;

for ii = 1:NumSimu
    
    if ii == 200000*simu_count
        clock_new = cputime;
        disp(['Simulation cycle is: ',num2str(ii),' & Global count = ',num2str(global_count),' & Required time is: ',num2str(clock_new-clock_old)]);
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
            Count_EventW = Count_EventW + 1;
        
        elseif ws_pick == 2 %chose magnitude and RI from section rupture list
   
        %for section ruptures, chose weighting between faults and
        %multifaults consistent with weight_SWM 
                tmp=sum(weight_SWM(1,1)+weight_SWM(1,3));
                ran_s = tmp*rand(1,1);   
                    if  ran_s < weight_SWM(1,1)
                        recurrence_para = recurrence_para2; %Chose parameters for unsegmented faults that are not multi faults
                        Count_EventS_1 = Count_EventS_1 + 1;
                    else
                        recurrence_para = recurrence_para3; %Chose parameters for unsegmented faults that are multi faults
                        Count_EventS_2 = Count_EventS_2 + 1;
                    end
                    
        elseif ws_pick ==3 %chose magnitude and RI from multifault rupture list
            
        %for multifault ruptures, chose weighting between faults and
        %sections consistent with weight_SWM 
                tmp=sum(weight_SWM(1,1)+weight_SWM(1,2));
                ran_m = tmp*rand(1,1); 
                    if  ran_m < weight_SWM(1,1)
                         recurrence_para = recurrence_para4; %Chose parameters for non mulitfaults that are faults 
                         Count_EventM_1 = Count_EventM_1 + 1;
                    else
                        recurrence_para = recurrence_para5; %Chose parameters for non mulitfaults that are also segmented
                        Count_EventM_2 = Count_EventM_2 + 1;
                    end
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
            
            % Top depth 
            % Assinged depths of 1, 3, or 5 km
            ran = rand(1,1);
            for jj = 1:num_TD
                if ran > sum(weight_TD(1,1:jj-1)) && ran <= sum(weight_TD(1,1:jj))
                    td = weight_TD(2,jj);
                    break
                end
            end
            
           rupture_id     = recurrence_para(num,3);  %Identify rupture id 
           global_count = global_count + 1;
                                                                      
            % Synthetic earthquake catalog
            % 1 ) Overall cycle
            % 2 ) Simulation cycle
            % 3 ) Occurrence time (since last event)
            % 4 ) Earthquake occurrence time (cumulative)
            % 5 ) Earthquake magnitude
            % 6 ) Top depth (km)
            % 7 ) Source ID 
            % 8 ) Section, whole fault, or multi fault rupture (1 to 3)
            EQCAT(global_count,1:8) = [global_count ii t_sum (ii-1)*t_limit+t_sum m td rupture_id ws_pick];
            
        end
        
    end
    
    Num_event(ii,1) = global_count - count_pre; %Num_event for each year in catalog
    
end

disp(['Total number of earthquakes: ',num2str(global_count)]);

%% Save to eqcat_test.mat. This is linked to EQCAT for psha by excel input file

save(char(syncat_name),'EQCAT','global_count','Num_event','recurrence_para1','recurrence_para2',...
    'recurrence_para3','recurrence_para4','recurrence_para5', 'num_sec1', 'num_multi_fault');


    
    
