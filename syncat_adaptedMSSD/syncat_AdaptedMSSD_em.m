%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%   Adapted Stochastic Event Catalog for all faults in Malawi   %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Same as Adatated Catalog but generates seperate catalog
% for each width and mf shape considered

% All fault source info from num_MSSD, where information for multi-faults
% and faults that are not part of multifault systems

clear all

mydir  = pwd; idcs   = strfind(mydir,'/');
addpath(mydir(1:idcs(end)-1));
addpath([mydir(1:idcs(end)-1) '/syncat_bg']);
addpath([mydir(1:idcs(end)-1) '/misc_functions']);

load MSSD_sources; load SZ_final;
load('syncat_PSHA_MSSD_input','bg_mmin','NumSimu','t_limit');

%% Select MSSD Fault info and magnitude recurrence characteristics

Mx = num_MSSD(:,12);%Max mag from Leonard scaling relationships, but this is varied
rigidity = 33*10^9; % N/m^2, set so consistent with Leonard (2010)
Mmin = bg_mmin;

%Parameter from Y&C 85 that governs magnitude increment over which fault
%hosts characteristic behaviour
DeltaM1 = zeros(length(Mx),1);

%Parameter from Y&C 85 that governs magnitude increment between G-R and
%characteristic behaviour
DeltaM2 = zeros(length(Mx),1); 

for i=1:length(Mx)
    if Mx(i)>6.5
        DeltaM1(i) = 1.0; DeltaM2(i) = 0.5; 
    else
        DeltaM1(i) = 0.5; DeltaM2(i) = 0.25;
    end
end
%Magnitude range in interval between Mmin and Mmax in G-R relation
% Varies depending on Mmax of each source
dm = 0.005;
mag_range_GR={};
for ff=1:length(num_MSSD)
mag_range_GR{ff} = ((Mmin-dm/2):dm:(Mx(ff)+0.15+dm/2))';
end

%Extract slip rate from num_MSSD in m/yr 
Slip_rate_option = num_MSSD(:,6)*10^-3;
    
%b-value for Malawi from GEM
Bvalue_option = [SZ_final(9).b+0.1 SZ_final(9).b SZ_final(9).b-0.1];
Bvalue_weight = [0.16 0.68 0.16];

Mmax_shift_option = [-0.15 0 0.15];  % Mmax shift 
Mmax_shift_weight = [0.3 0.6 0.1];

%% Cases for recurrence relationships
% Explores all combinations of slip rates, b value, and m_max shift, and
% the weighting for each of these combinations
count_case = 0;
RecurrenceVar = {};

for ff = 1:length(num_MSSD)
%for ii = 1:length(Slip_rate_option)
    for jj = 1:length(Bvalue_option)
        for kk = 1:length(Mmax_shift_option)
            count_case = count_case + 1;
            %RecurrenceVar(count_case,1:5) = [count_case Slip_rate_option(ii) Bvalue_option(jj) Mmax_shift_option(kk) Slip_rate_weight(ii)*Bvalue_weight(jj)*Mmax_shift_weight(kk)];
            RecurrenceVar{ff}(count_case,1:5) = [count_case Slip_rate_option(ff) Bvalue_option(jj) Mmax_shift_option(kk) Mmax_shift_weight(kk)*Bvalue_weight(jj)];
            %RecurrenceVar{ff}(count_case,1:4) = [count_case Slip_rate_option(ff) Bvalue_option(jj) Bvalue_weight(jj)];
        %end
        end
    end
   count_case = 0;
end

%% Run fault recurrence model for multiple faults and different fault widths
%Mmax still not considered an uncertainity

%Save variables in cell arrays for each fault
alpha_var1={}; fm_char_var1={}; fm_exp_var1={};%variables for width=1
alpha_var2={}; fm_char_var2={}; fm_exp_var2={};%variables for width=2
YC85SeismicMoment_exp={}; YC85SeismicMoment_char={}; %variables to store moment rates in

for ff= 1:length(num_MSSD) %loop through each fault/multi_fault

    SeismicMoment_char = zeros(height(RecurrenceVar{ff}),3);
    SeismicMoment_exp  = zeros(height(RecurrenceVar{ff}),3);

    weighted_GR_char = zeros(length(mag_range_GR{ff}),2);
    weighted_GR_exp  = zeros(length(mag_range_GR{ff}),2);   
  
    for ii = 1:2
    
        if ii == 1
            Af = num_MSSD(ff,9)*10^6;%Width controlled by source length
        elseif ii == 2
            Af = num_MSSD(ff,11)*10^6;%Width controlled by 35 km thick seismogenic layer
        end

        for jj = 1:height(RecurrenceVar{ff})
        
        %Interpreted differently from Covertino as mc = mmax
        %In other words magnitude range = MX-delta 2, not M+delta 2  
        Mmax = Mx(ff) + RecurrenceVar{ff}(jj,4);%vary mmax +/- 0.15
    
        % Input variables of para_GR:
        % 1) Rigidity (Nm)
        % 2) Fault plane area Af (m^2) for fault ff and case jj
        % 3) Slip rate S (m/year) for fault ff and case jj
        % 4) Gutenberg-Richter slope parameter b
        % 5) Maximum magnitude Mmax for fault ff
        % 6) Minimum magnitude Mmin
        % 7) Magnitude range for exponential distribution part DeltaM1
        % 8) Magnitude range for characteristic distribution part DeltaM2
        
        para_GR = [rigidity Af RecurrenceVar{ff}(jj,2) RecurrenceVar{ff}(jj,3) Mmax Mmin DeltaM1(ff) DeltaM2(ff)];
        
        %Run YC1985 function where inputs are: 1) para_GR that includes SR and Mmax info,
        %2) mag range for GR relation, and 3) tick labels for a G-R plot
        
        %Outputs are activity rate and pdf and cdf for GR and
        %characteristic behaviours
        [alpha_var(jj,1:3),fm_char_var(1:length(mag_range_GR{ff}),3*(jj-1)+1:3*jj),fm_exp_var(1:length(mag_range_GR{ff}),3*(jj-1)+1:3*jj)] = ...
            characteristic_magnitude_YC1985(para_GR,mag_range_GR{ff},[],[6 8 10^-5 10^-3]);
        
        %weighted mean of frequency for each magnitude for case where width =1 and width =2
        
        %Note weighted average cuts out at mmax-mshift (i,e lower bound of
        %Mmax). Not true MFD can't be used for moment rate calc
          weighted_GR_char(:,ii) = weighted_GR_char(:,ii) + RecurrenceVar{ff}(jj,5)*fm_char_var(:,3*jj);
          weighted_GR_exp(:,ii)  = weighted_GR_exp(:,ii)  + RecurrenceVar{ff}(jj,5)*fm_exp_var(:,3*jj);
          
        for kk = 1:length(mag_range_GR{ff})
            if fm_char_var(kk,3*jj) > 0
                %moment rate for each case from assessing rate from each
                %magnitude interval
                SeismicMoment_char(jj,ii) = SeismicMoment_char(jj,ii) + (fm_char_var(kk,3*jj)-fm_char_var(kk+1,3*jj))*(10^(1.5*mag_range_GR{ff}(kk)+9.05));
 
            end
            if fm_exp_var(kk,3*jj) > 0
                SeismicMoment_exp(jj,ii) = SeismicMoment_exp(jj,ii) + (fm_exp_var(kk,3*jj)-fm_exp_var(kk+1,3*jj))*(10^(1.5*mag_range_GR{ff}(kk)+9.05));
            end
            
        end
        
        YC85SeismicMoment_exp{ff}(jj,ii)=SeismicMoment_exp(jj,ii);
        YC85SeismicMoment_char{ff}(jj,ii)=SeismicMoment_char(jj,ii);
        
        end  %end jj loop for different b-value and slip rate variations
  
    
    % Assign fault ID's to data
    alpha_var = horzcat(alpha_var,ones(height(alpha_var),1)*num_MSSD(ff,1));
    
    % Save the fault information for PSHA
    if ii == 1 %Variables for fault limited width
        alpha_var1{ff} = alpha_var; fm_char_var1{ff} = fm_char_var; fm_exp_var1{ff} = fm_exp_var;
        
    elseif ii == 2 %Variables for layer-limited width
        alpha_var2{ff} = alpha_var; fm_char_var2{ff} = fm_char_var; fm_exp_var2{ff} = fm_exp_var;
        
    end
    
   clear alpha_var fm_char_var fm_exp_var 
 

end %end ii loop for width
    clear SeismicMoment_exp  SeismicMoment_char
end % end loop for fault

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%   SYNTHETIC CATALOG  %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Combine cell arrays from each fault for use in syncat
% Last column gives fault id
alpha_var1_comb = cell2mat(alpha_var1');
alpha_var2_comb = cell2mat(alpha_var2');
%Based on codes for PSHA LRVF_DMF PSHA
%However, in this case just using it to make synthetic event catalog

%Set-up

rnd_gen_index = 0;
if rnd_gen_index == -2
    rand('state',sum(clock*100));
else
    rand('state',rnd_gen_index);
end

disp(['Total duration of simulation: ',num2str(NumSimu*t_limit),' years']);

% Layer 3 : B value shift and Mmax shift
recurrence_para_weight = RecurrenceVar{1}(:,5)';

fault_plane_visual_check = [];
% fault_plane_visual_check = 909

clock_old = cputime;

count_case=0;
syncat_opt=zeros(1,2);

%syncat opt: column 1 = fault width case
%syncat opt: column 2 = char (1) or GR (2)

for ww=1:2 %select of mfd
    for mm =1:2 
        count_case=count_case+1;
        syncat_opt(count_case,1:2)=[ww mm];
        
        %final length of catalog unknown, but can assume it will be less than 4x
        %length of simulation. Set length now to save time, and remove zeros later
        StochasticEventCatalog{count_case}=zeros(4*NumSimu,8);
    end
end

%% Run simulations

for ss=1:length(syncat_opt)

    global_count = 0; Num_event_adapted_tmp=zeros(NumSimu,1);
    tmp_count = 1;
    
    fault_width=syncat_opt(ss,1); %select fault width case
    recurrence_type=syncat_opt(ss,2); %select recurrence_type
    
    %chose activity rates and MFD based on fault width case   
    if fault_width==1
        alpha_var_comb=alpha_var1_comb;
        fm_char_var=fm_char_var1;
        fm_exp_var=fm_exp_var1;
    else
        alpha_var_comb=alpha_var2_comb;
        fm_char_var=fm_char_var2;
        fm_exp_var=fm_exp_var2;
    end
    
for ii = 1:NumSimu %Run event catalog for duration set by N 
    
    if ii == tmp_count*50000
        tmp_count = tmp_count + 1;
        clock_new = cputime;
        disp(['Current simulation cycle : ',num2str(ii),' out of ',num2str(NumSimu),' & Required time is: ',...
            num2str(clock_new-clock_old), '. Catalog is ', num2str(ss), ' out of ', num2str(length(syncat_opt))]);
        clock_old = clock_new;
    end
  
    % Logic tree sampling
    ran_LT = rand(1,1);%1 random numbers between 0 and 1  
   
    % Recurrence parameters (Bvalue, and Mmax shift)
    for jj = 1:length(recurrence_para_weight)
        if ran_LT <= sum(recurrence_para_weight(1:jj))
            recurrence_para = jj; 
            break
        end
    end     
       
    % Characteristic events
    if recurrence_type == 1
        %activity rate for characteristic events, includes both 
        %exponential and non-exponential part of MFD
        %index from alpha_var variable at intervals of number of faults
        alpha_info_tmp = alpha_var_comb((recurrence_para:9:length(alpha_var_comb)),1:2);
        alpha_info = sum(alpha_info_tmp,2);
            
    clear alpha_info_tmp
            
    %G-R events    
    elseif recurrence_type == 2
            
        %Chose activity rate and reccurence rate from exp MDF
        %Only one activty rate to consider
        alpha_info  = alpha_var_comb((recurrence_para:9:length(alpha_var1_comb)),3); 
    end
     
        %Occurrence times estimated from 1/alpha rate where alpha rate is
        %the fault's activity rate. In other words, the frequency of events
        %over Mmin (4.5), so 1/alpha is the recurrence interval of events
        %>Mmin
        
    for ff=1:length(alpha_info)
        %need inverse of activity rate (ie recurrence interval) to
        %set as mean in expinv
        alpha_info(ff,2)=1/alpha_info(ff,1);
    end
        
    occurrence_time = zeros(1,length(alpha_info));
    %Uses expinv, where event is taken if occurrence time is below 1
    %Needs to be done for each fault
    occurrence_time = expinv(rand(1,length(occurrence_time)),alpha_info(:,2)'); 
        
    if min(occurrence_time) <= t_limit
            
       f_indx = find(occurrence_time <=t_limit);%index faults for which event occurred
       global_count = global_count + length(f_indx);
            
       for rr=1:length(f_indx) %for number of faults with events in this simulation cycle
           if recurrence_type == 1 %if char MFD
              %index mag from fault specific GR MFD CFD
              recurrence_info1 = fm_char_var{f_indx(rr)}(:,3*(recurrence_para-1)+2);
              m{rr} = mag_range_GR{f_indx(rr)}(find(rand(1,1)>recurrence_info1,1,'last')+1) + rand(1,1)*dm - dm/2;
           elseif recurrence_type == 2 %if exp MFD
              %index mag from fault specific char MFD CFD
              recurrence_info1 = fm_exp_var{f_indx(rr)}(:,3*(recurrence_para-1)+2);
              m{rr} = mag_range_GR{f_indx(rr)}(find(rand(1,1)>recurrence_info1,1,'last')+1) + rand(1,1)*dm - dm/2;
           end %end if loop
       
    % StochasticEventCatalog
    % 1    ) Overall cycle
    % 2    ) Simulation cycle
    % 3    ) Occurrence time
    % 4    ) Earthquake magnitude
    % 5    ) Fault model (1 = partial width versus 2 = full width)
    % 6    ) MSSD ID
    % 7    ) Recurrence type (1 = characteristic versus 2 = exponential)
    % 8    ) Reccurrece parameter
     
         StochasticEventCatalog{ss}(global_count-length(f_indx)+rr,1:8) = [global_count-length(f_indx)+rr ii occurrence_time(f_indx(rr)) m{rr} fault_width num_MSSD(f_indx(rr),1) recurrence_type recurrence_para];
       end %end rr loop
            
   Num_event_adapted_tmp(ii,1) = length(f_indx); %Num_event for each year in catalog 
          
    else
    %Record if no events during that simulation
    
   Num_event_adapted_tmp(ii,1) = 0;
   
    end
  
   
   clear m
            
end %end simulation cycle

    clear alpha_var_comb fm_char_var fm_exp_var recurrence_type_fault_width 
    
    %Remove unnecessary zeros;
    StochasticEventCatalog{ss}=StochasticEventCatalog{ss}(1:sum(Num_event_adapted_tmp),:);
    Num_event_adapted{ss}=Num_event_adapted_tmp;
    
end %end ss loop
%% Save or load file for analysis

save('MSSD_Catalog_Adapted_em','StochasticEventCatalog','NumSimu','Num_event_adapted','t_limit','YC85SeismicMoment_char','YC85SeismicMoment_exp',...
    'RecurrenceVar','dm','mag_range_GR','fm_char_var1','fm_char_var2',...
    'fm_exp_var1','fm_exp_var2','syncat_opt')
    