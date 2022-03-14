%% SYNCAT ANALYSIS

%Analyses and tests EQCAT generated from syncat_MSSD

load eqcat_test %Load direct MSSD fault record
mydir  = pwd; idcs   = strfind(mydir,'/');
addpath(mydir(1:idcs(end)-1));

load MSSD_sources
load('syncat_PSHA_MSSD_input','NumSimu','weight_SWM','t_limit')

%% Load syncat_bg for comparison
addpath([mydir(1:idcs(end)-1) '/syncat_bg']);
load syncat_bg

%Remove events outsie Malawi and those in Malawi Mw>7, as these events are 
%can instead be associated with faults
for i=1:length(bg_catalog)
     if bg_catalog(i,1) > 0 %positive distances equal outside PSHA region
        bg_catalog(i,5) = 0;
     else
        if bg_catalog (i,1)<0 && bg_catalog(i,4) >7.0 
           bg_catalog(i,5) = 0;
        end
     end
end

bg_catalog(bg_catalog(:,5)==0,:)=[];

%% Check the MFD of the simulated earthquakes and compare with syncat_bg

mag_check = 4.5:0.05:8.5;
Rate_bg=zeros(length(mag_check),2);  Rate=zeros(length(mag_check),1); Rate_c = zeros(length(mag_check),1);


for ii = 1:length(mag_check)
      %Cumulative annual number of events greater than each mag interval (i)
      Rate_bg(ii,1) = length(find(bg_catalog(:,4)>=mag_check(ii)))/(NumSimu*t_limit);
      %Annual number of events for discrete mag interval(i)
      Rate_bg(ii,2) = length(find(bg_catalog(:,4)>=mag_check(ii) & bg_catalog(:,4)<mag_check(ii)+0.1))/(NumSimu*t_limit); 
      Rate(ii,1) = length(find(EQCAT(:,5)>=mag_check(ii)))/(NumSimu*t_limit); %Bin events cumulatively by magnitude at intervals of 0.1 
      Rate_c(ii,1) = Rate(ii,1)+Rate_bg(ii,1); %Combined records
end
%Plot freq-mag of EQ cat with instrumental G-R    

figure(1);
semilogy(mag_check,Rate(:,1),'LineWidth',1.5); hold on;
semilogy(mag_check,Rate_bg(:,1),'r-','LineWidth',1.5); hold on;
semilogy(mag_check,Rate_c(:,1),'k-','LineWidth',1.5); hold on;

xlabel('Magnitude'); ylabel('Annual frequency of exceedance');
axis square; grid on; axis([5 9 10^-5 10^0]);
legend('MSSD sources','Bg Catalog','Combined');
set(gca,'FontSize',13.5);


%% Seismic Moment Rate Calcs

%EQCAT
for i =1:length(EQCAT)
    s_Mo(i)=10^(EQCAT(i,5)*1.5+9.05); %Calcuate seismic moment of each event
end


TotalMo = sum(s_Mo); %Sum seismic moment of all events from each catalog
Moment_Rate=TotalMo/NumSimu; %Divide total seismic moment by duration of record

%bg_catalog
for i =1:length(bg_catalog)
    s_Mo_bg(i)=10^(bg_catalog(i,4)*1.5+9.05); %Calcuate seismic moment of each even
end

TotalMo_bg = sum(s_Mo_bg); %Sum seismic moment of all events from each catalog
Moment_Rate_bg=TotalMo_bg/NumSimu; %Divide total seismic moment by duration of record

Mo_c = Moment_Rate_bg+Moment_Rate;

%% Check moment rate of faults in syncat with rate from database

syncat_fault_mr=zeros(length(num_fault),1);
syncat_sec_mr=zeros(length(num_sec1),1);
syncat_multi_fault_mr=zeros(height(num_multi_fault),1);
count=0;
for i=1:length(num_fault)
    
    if isempty(find(num_fault(i,1)==recurrence_para2(:,3))) ==0 && isempty(find(num_fault(i,1)==recurrence_para4(:,3))) ==0
        %Find unsegmented faults that is not part of multi-fault system
        %Moment rate from dividing by the full length of the simulation
        syncat_fault_mr(i)=sum(10.^(9.05+1.5*(EQCAT(find(num_fault(i,1)==EQCAT(:,7)),5))))/NumSimu;
        
    elseif isempty(find(num_fault(i,1)==recurrence_para2(:,3))) ==1 && isempty(find(num_fault(i,1)==recurrence_para3(:,3))) ==1 ...
        && isempty(find(num_fault(i,1)==recurrence_para4(:,3))) ==1
        %Find segmented faults that is part of multi-fault system
        syncat_fault_mr(i)=sum(10.^(9.05+1.5*(EQCAT(find(num_fault(i,1)==EQCAT(:,7)),5))))/(NumSimu*weight_SWM(1,1));      
        
    elseif isempty(find(num_fault(i,1)==recurrence_para2(:,3))) ==1 && isempty(find(num_fault(i,1)==recurrence_para4(:,3))) ==0 ...
        && isempty(find(num_fault(i,1)==recurrence_para5(:,3))) ==1
        %Find segmented faults that is not part of multi-fault system
        syncat_fault_mr(i)=sum(10.^(9.05+1.5*(EQCAT(find(num_fault(i,1)==EQCAT(:,7)),5))))/(NumSimu*(weight_SWM(1,1)+weight_SWM(1,3)*weight_SWM(1,1)));   
   
    elseif isempty(find(num_fault(i,1)==recurrence_para2(:,3))) ==0 && isempty(find(num_fault(i,1)==recurrence_para3(:,3))) ==1 ...
        && isempty(find(num_fault(i,1)==recurrence_para4(:,3))) ==1
        %Find unsegmented faults that is part of multi-fault system
        syncat_fault_mr(i)=sum(10.^(9.05+1.5*(EQCAT(find(num_fault(i,1)==EQCAT(:,7)),5))))/(NumSimu*(weight_SWM(1,1)+weight_SWM(1,2)*(weight_SWM(1,1)/(weight_SWM(1,1)+weight_SWM(1,3)))));   
    end
end


for i=1:length(num_sec1)
     if  num_sec1(i,18)>0
         %Find if segment is part of multi fault system
        syncat_sec_mr(i)=sum(10.^(9.05+1.5*(EQCAT(find(num_sec1(i,1)==EQCAT(:,7)),5))))/(NumSimu*weight_SWM(1,2));
     else
         %If not:
        syncat_sec_mr(i)=sum(10.^(9.05+1.5*(EQCAT(find(num_sec1(i,1)==EQCAT(:,7)),5))))/(NumSimu*(weight_SWM(1,2)+weight_SWM(1,3)*(weight_SWM(1,2)/(weight_SWM(1,1)+weight_SWM(1,2)))));
     end
end

for i=1:height(num_multi_fault)
     if  isempty(find(num_multi_fault(i,1)== num_sec1(:,18)))==0
         %Find if multi_fault part of segment
        syncat_multi_fault_mr(i)=sum(10.^(9.05+1.5*(EQCAT(find(num_multi_fault(i,1)==EQCAT(:,7)),5))))/(NumSimu*weight_SWM(1,3));
     else
         %If not:
        syncat_multi_fault_mr(i)=sum(10.^(9.05+1.5*(EQCAT(find(num_multi_fault(i,1)==EQCAT(:,7)),5))))/(NumSimu*(weight_SWM(1,3)+weight_SWM(1,2)*(weight_SWM(1,3)/(weight_SWM(1,1)+weight_SWM(1,3)))));
     end
end

% Combine moment rate from sections that are identical except in their participation in different fault rupture

%LIST OF DIFFERENT SECTIONS WITH SAME MAGNITUDES
%Liwawadzi South & Sani Central 
%Malombe South & NB Fault 9a&b
%Mwanza Thombani, Wamkurumadzi South-1 & Chlingali South
%Sani North, Mbiri-1 Central,Kaporo-1 North
%BMF  Matakataka-1&2, Lweya Central
%Ruo North & Metangula Lake 1&2
%Zomba Chingale Stream & Phirilanyama-2
%Thyolo Mbewe & Metangula North 2a
%Central Basin Fault 20 North, Wovwe-2 North, St Mary South
%Bwanje North & South Basin Fault 8 South
%Chirobwe Ncheu-2 Dzonze South, Malombe Chimwalira, South Basin Fault 7b South
%Chingale Step Namitembo, BMF Kasinje 1&2, BMF Lithipe-2
%South Basin Fault 13c South & Central Basin Fault 19 North
%Mwanza Condedezi & South Basin Fault 7 North a&b
%Zomba North & Makanjira Malindi
%South Basin Fault 13b South & Liwawazi North
%Chirobwe Ncheu Fault Livuzeli 1&2 & Malombe Mpale
%Mbiri-1 North & Lipichilli-3
%Thyolo North 1&2 & South Basin Fault 7a

match_sec=[94,24,79,3,44,98,70,95,149,128,8,17,96,91,160,134,140,162,49,123,58,78,121,...
    22,73,117,133,4,19,85,25,115,76,138,143,119];

for i=1:length(num_sec1)
    %find secs with matching magnitudes
    
   if length(find(num_sec1(i,21)==num_sec1(:,21)))>=2  && isempty(find(i==match_sec))==1
       
        tmp_indx=find(num_sec1(i,21)==num_sec1(:,21)); 
       
        %Correction for lengths of:
        %NB-F 9a and 9b
        %Metangula Lake 1&2
        %South Basin Fault 7 North a&b
        if isempty(find(i==[167,165,87,90,120,122]))==0
        tmp_indx(1)=[];
        %Correction for lengths:
        %of Mbiri with Kaporo & Sani
        elseif isempty(find(i==[141,144,62,70]))==0
        tmp_indx([1,4])=[];
        %Correction for lengths:
        %BMF  Matakataka-1&2, Lweya Central
        %Thyolo Mbewe & Metangula North 2a
        %Chirobwe Ncheu Fault Livuzeli 1&2 & Malombe Mpale
        %Thyolo North 1&2 & South Basin Fault 7a
        elseif isempty(find(i==[60,68,10,14,51,55,9,13]))==0
        tmp_indx(3)=[];
        end
        
        for ii=1:length(tmp_indx)
         
        num_sec1(tmp_indx(1),26)=sum(num_sec1(tmp_indx,26));
        syncat_sec_mr(tmp_indx(1))=sum(syncat_sec_mr(tmp_indx)); 
            if ii>1 
            num_sec1(tmp_indx(ii),26)=0; syncat_sec_mr(tmp_indx(ii))=0;
            end
        end
        clear tmp_indx
   end
      
end

% combine moment rate from faults that are identical except in
% their participation in different multifault rupture

%LIST OF DIFFERENT FAULTS WITH SAME MAGNITUDES
%Bwanwie (25) & Lipichilli-3 (87)
%South Fault Basin 6 (58) & Central Basin Fault 7 (81)
%Mbiri-2 (98) & Usisya-North (70)

match_fault=[25,87,58,81,98,70];

for i=1:length(num_fault)
    %match these faults based on equal mag estimates
    %Note 
    %length, so are exluded from this comparison
    
   if (length(find(num_fault(i,18)==num_fault(:,18)))>=2) && isempty(find(i==match_fault))==1
       
       tmp_indx=find(num_fault(i,18)==num_fault(:,18));
       %If fault has been corrected before or has an unequal chance of being part multifault system
       %don't include.
        for ii=1:length(tmp_indx)
         
        num_fault(tmp_indx(1),25)=sum(num_fault(tmp_indx,25));
        syncat_fault_mr(tmp_indx(1))=sum(syncat_fault_mr(tmp_indx)); 
            if ii>1
            num_fault(tmp_indx(ii),25)=0; syncat_fault_mr(tmp_indx(ii))=0;
            end
        end
        clear tmp_indx
        
   end
     
end

%Remove duplicate sections and faults
num_sec1(find(num_sec1(:,26)==0),:)=[]; num_fault(find(num_fault(:,25)==0),:)=[];
syncat_sec_mr(find(syncat_sec_mr==0))=[]; syncat_fault_mr(find(syncat_fault_mr==0))=[]; 


%% PLOT FIGURE AND SAVE CALCULATIONS FOR ANALYSIS IN 'MSSD_COMB'

figure(2);

plot([13.5 17.5],[13.5 17.5],'k--',... 
log10(num_sec1(:,26)),log10(syncat_sec_mr),'rx',...
log10(num_fault(:,25)),log10(syncat_fault_mr),'bx',...
log10(num_multi_fault(:,15)),log10(syncat_multi_fault_mr),'gx','MarkerSize',12,'LineWidth',1.5);
axis([13.5 17.5 13.5 17.5]); axis square; hold on;
xlabel('MSSD Moment Rate (log Nm/yr)'),ylabel(['Stochastic Event Catalog' newline 'Moment Rate (log Nm/yr)']); hold on;
legend('y=x',['section ruptures n = ' num2str(length(num_sec1))],['fault ruptures n=' num2str(length(num_fault))],...
    ['multi-fault ruptures n =' num2str(length(num_multi_fault))] ,'Location','southeast');
grid on; set(gca,'fontsize',13)

save('syncat_comparison','syncat_sec_mr','syncat_fault_mr','syncat_multi_fault_mr','num_sec1','num_fault','num_multi_fault');

