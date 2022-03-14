%                                              %
% M-function for obtaining annual maximum data %
%                                              %
    
%Parameters are set in PSHA_MSSD, and are vectors of 0, where the number
%of rows is equal to the number of events within each group

% Iteratively takes ground motions from previous groups and from
% most recent siumations to find highest ground motions

% Ground motion
% EPS
% Distance between fault and site (4 possible choices)
% Mag of event
% TIME (1) num of event, (2) absolute time of event
% GMPE ID
% SZ (length 3): Fault ID, catalog pick (segmented, fault, multi_fault)


function [AMAX_GM,num_taken] = annual_max_MSSD(GM,EPS,DIST,MAG,TIME,GMPEID,SZ,num_T,group_count,num_need,num_taken,AMAX_GM)

n_bin = 1;

% PUTS GM RECORDS FROM EQCAT INTO VECTOR WITH LENGTH OF GROUP_SIZE (I.E 10000) SORTED BY VIBRATION PERIOD %
%{
for jj = 1:num_T
    max_GM(1,1:10,jj) = [GM(jj,1) EPS(jj,1) DIST(4,1) MAG(1) TIME(1,1) TIME(2,1) GMPEID(1) DIST(2,1) SZ(1,1) SZ(2,1)];
end

%}

max_GM=zeros(length(MAG),10,num_T);

for kk = 1:length(MAG) %loop through num_talem
    
    %{
    Option to find events that occurred in the same simulation cycle (i.e. year), and take only
    the largest event. Option removed folloing discussion with Katsu (see email 19/03/21)
    if TIME(2,kk) == max_GM(n_bin,5,1)% && floor(TIME(3,kk)) == floor(max_GM(n_bin,6,1)) %For event times that match 
        
        for jj = 1:num_T
            if max_GM(n_bin,1,jj) < GM(jj,kk) %For simulations in syncat when two events, selects event with largest GM
                max_GM(n_bin,1:10,jj) = [GM(jj,kk) EPS(jj,kk) DIST(4,kk) MAG(kk) TIME(1,kk) TIME(2,kk) GMPEID(kk) DIST(2,kk) SZ(1,kk) SZ(2,kk)];
                kk_tmp(kk)=kk;
            end
        end
        
        else
     %}
        n_bin = n_bin + 1;
        
        for jj = 1:num_T
            max_GM(n_bin,1:10,jj) = [GM(jj,kk) EPS(jj,kk) DIST(4,kk) MAG(kk) TIME(1,kk) TIME(2,kk) GMPEID(kk) DIST(2,kk) SZ(1,kk) SZ(2,kk)];
        end           
    %end 
end
          
% Store only necessary samples for further analyses
% takes whatever is smaller, all events (num_taken) or num_need
num_taken = num_taken + n_bin;
if group_count == 1
    
    AMAX_GM = max_GM;
    
else
    
    tmp_AMAX_GM = zeros(min(num_need,num_taken),10,num_T);
    for j = 1:num_T
        %combine with AMAX_GM from previous groups
        tmp(:,:) = [AMAX_GM(:,:,j); max_GM(1:n_bin,:,j)];
        tmp = -sortrows(-tmp,1); %Sort rows in descending order of GM
        %only take highest ground motions from most recent simulations (max_GM)
        %and previous runs (AMAX_GM)
        tmp_AMAX_GM(1:min(num_need,num_taken),1:10,j) = tmp(1:min(num_need,num_taken),1:10);
    end
    
    clear AMAX_GM
    AMAX_GM = tmp_AMAX_GM;
    clear tmp tmp_AMAX_GM
    
    num_taken = min(num_taken,num_need);
    
end

% Returns matrix with rows for given spectral acceleration and event:

% 1 Ground Motion
% 2 EPS
% 3 Distance (Rhypo)
% 4 Event magnitude
% 5 Time(1) (event num)
% 6 Time(2) (simulation number)
% 7 GMPE index
% 8 Distance (Rjb)
% 9 Source type (sec, fault, mf)
% 10 Source id


