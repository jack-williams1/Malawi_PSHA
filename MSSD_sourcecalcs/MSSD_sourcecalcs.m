%% CALCULATE FAULT SLIP RATES AND EARTHQUAKE RECURRENCE INTERVALS 
% IN THE MSSD THROUGH MONTE CARLO SIMULATION THROUGH LOGIC TREE
% SIMULATION RESULTS CAN ALSO BE USED TO DERIVE EPESTEMIC UNCERTAINITY

%load and sort data

%load source inputs

[num_raw,txt_raw,raw_raw]=xlsread('MSSD.xlsx',1);

%load hwf correction

[num_hwf,txt_hwf,raw_hwf]=xlsread('HangingWallFlexureInputs.xlsx','A5:AE11');

%load north basin fault data to use instead
addpath('nbasin_faultdata'); nbasin_faults_sr=readtable('nb_slipratetable_id.csv');
nbasin_faults=table2array(nbasin_faults_sr);

mydir  = pwd; idcs   = strfind(mydir,'/');
addpath(mydir(1:idcs(end)-1));
load MSSD_sources

%load logic tree weightings
load('syncat_PSHA_MSSD_input','dip_w','hwf_w','ext_w')

%Remove linking sections
%i.e those sections who don't rupture on their own and have a moment rate set as '0'
tmp_bb=find(num_sec(:,26)); num_sec1=num_sec(tmp_bb,:); clear tmp_bb

%link rows in num_raw with MSSD_sources variables and basin

MSSD_id = vertcat(num_fault(:,1),num_sec1(:,1),num_multi_fault(:,1));
basin_indx=zeros(length(MSSD_id),1);

for i=1:length(MSSD_id)
    
    if MSSD_id(i,1) >600
        MSSD_id(i,2)=find(MSSD_id(i,1)==num_raw(:,83),1);
    elseif MSSD_id(i,1) >300
        MSSD_id(i,2)=find(MSSD_id(i,1)==num_raw(:,81),1);
    else
        MSSD_id(i,2)=find(MSSD_id(i,1)==num_raw(:,82),1);
    end
    
    %basin index for hanging wall flexure correction
    if isempty(find(string(txt_raw(MSSD_id(i,2)+2,3))==string(txt_hwf(:,1))))==0    
        basin_indx(i)=find(string(txt_raw(MSSD_id(i,2)+2,3))==string(txt_hwf(:,1)));
    end
end

%% Logic tree set-up
% Layer 1 : Fault Dip
% Layer 2 : Rift extension rate
% Layer 3 : Azimuth of rift extension
% Layer 4 : HW flexure correction
% Layer 5 : Distribution of rift extension between border and intrabasin fault
% Layer 6 : C1 parameter of Leonard (2010)
% Layer 7 : C2 parameter of Leonard (2010)

% Weightings are set in PSHA_MSSD_input

dip = [num_raw(:,13) num_raw(:,14) num_raw(:,15)];

ext_r = [num_raw(:,32) num_raw(:,33)];
ext_az = [num_raw(:,35) num_raw(:,36)];
hwf_f = [num_hwf(:,28) num_hwf(:,29) num_hwf(:,30)];
ext_d = [num_raw(:,29) num_raw(:,30) num_raw(:,31)];
c1 = [12 17.5];
c2 = [1.5 3.8];

sec_source_az = num_raw(:,8); flt_source_az = num_raw(:,10); m_flt_source_az = num_raw(:,12);

%Update id's of lake border faults if new faults added
%Include Lweya and Kuvuzi as in footwall of Usisya Fault
%Include SB Fault 3, 6, 8, 15 as in Central-South basin transition zone
lb_flt_id = [105:1:116,129:1:142,147:1:157,190:1:192,...
            337,338,351:1:358,363,366:1:374,406,...
            608:1:610,612];
        
%% Run simulations
 
n_sim =10000; crust_mu= 33*10^9; % N/m^2, set so consistent with Leonard (2010)

mean_sr =zeros(length(MSSD_id),1); std_sr=zeros(length(MSSD_id),1);
mean_r=zeros(length(MSSD_id),1); std_r=zeros(length(MSSD_id),1);
range_sr =zeros(length(MSSD_id),2); range_r =zeros(length(MSSD_id),2);
constraint_test=zeros(length(MSSD_id),1);

simu_count = 1; nbindx=0;



for ss=1:length(MSSD_id)

    if ss== 40*simu_count
        disp(['MSSD source: ',num2str(ss),' out of ',num2str(length(MSSD_id))]);
        simu_count = simu_count + 1;
    end

sslip_rate =zeros(n_sim,1); 
sed =zeros(n_sim,1); ri =zeros(n_sim,1);
    
  id=MSSD_id(ss,1); ii=MSSD_id(ss,2); 

 if id>600 %Chose parameters from multi_fault
     source_az=m_flt_source_az(ii);
     len=num_raw(ii,11); area1 = num_raw(ii,25); area2 = num_raw(ii,26);
   elseif id>300 %Chose parameters from fault
     source_az=flt_source_az(ii);
     len=num_raw(ii,9); area1 = num_raw(ii,21); area2 = num_raw(ii,22);
   else %Chose parameters from section
     source_az=sec_source_az(ii);
     len=num_raw(ii,7); area1 = num_raw(ii,21); area2 = num_raw(ii,22);
     l_sec_indx=find(id==num_sec(:,1));
 end

 
 %Derive 10000 picks for continuous functions
 
 % Rift extension rate [continuous]
 ext_r_pick=normrnd(ext_r(ii,2),ext_r(ii,2)-ext_r(ii,1),n_sim,1);
 
 % Rift extension azimuth [continuous]
 ext_az_pick=normrnd(ext_az(ii,2),ext_az(ii,2)-ext_az(ii,1),n_sim,1);
 
 % C1 from Leonard 2010  [continuous]
 c1_pick = normrnd(c1(2),c1(2)-c1(1),n_sim,1);
    
 % C2 from Leonard 2010  [continuous, log normal] 
 c2_pick = (lognrnd(log(c2(2)),log(c2(2)-c2(1)),n_sim,1))/10^5;
 
 %For non-indexed southern basin, and border faults, hanging-wall flexure correction =1
 if basin_indx(ss)==0 || ismember(id,lb_flt_id)==1
     hwf_opt = [1 1 1];
 else 
     hwf_opt = hwf_f(basin_indx(ss),:);
 end
 
% Sample logic tree
for i=1:n_sim    
   
    %3 random numbers between 0 and 1 for each layer of logic tree
    ran_LT = rand(1,3);
    
    %index fault id in num_raw

  % Fault dip [discrete]
    for jj = 1:length(dip_w)
        if ran_LT(1) <= sum(dip_w(1:jj))
            dip_pick = jj; %If ran_LT(1)>0.9 jj=2
            break
        end
    end
    
    % Rift extension weight [discrete]
    for jj = 1:length(ext_w)
        if ran_LT(2) <= sum(ext_w(1:jj))
            ext_w_pick = jj; %If ran_LT(1)>0.9 jj=2
            break
        end
    end
   
    % Hanging wall flexure [discrete]
    % =1 f
    for jj = 1:length(hwf_w)
        if ran_LT(3) <= sum(hwf_w(1:jj))
            hwf_pick = jj; %If ran_LT(1)>0.9 jj=2
            break
        end
    end
    
    %Use picks to derive fault slip rate using eq. 3 in Williams et al (2021)
    sslip_rate(i) = ext_d(ii,ext_w_pick)*ext_r_pick(i)*hwf_opt(hwf_pick)*abs(cos(((source_az+90)...
    -ext_az_pick(i))*pi/180))/(cos(dip(ii,dip_pick)*pi/180));

    %Enforce postive slip rate
    if sslip_rate(i)<0
        sslip_rate(i)=0.001;
    end

    %IF FAULT SLIP RATE PREVIOUSLY MEASURED FROM OFFSET 75 KA REFLECTOR
    %USE THIS SLIP RATE INSTEAD
 
    if ismember(id,nbasin_faults(:,1))==1 || ... %index NB basin faults
     id<300 && ismember(num_sec(l_sec_indx,15),nbasin_faults(:,1))==1 || ...  %index NB basin fault secs
     id>600 && ismember(id,nbasin_faults(:,2))==1 ...  %index NB basin multifaults
    
     %if fault's first simulation, index in nbasin fault table
        if i==1
           nbindx=nbindx+1;
           slip_rate_nb.id{nbindx}=id;
           if id>600
           l_sr_indx=find(id==nbasin_faults(:,2),1);
           
           %find participating fault areas for multifaults that include
           %slip rates from faults that do and do not offset reflector
     
            if length(find(id==nbasin_faults(:,2)))==1
               flt_indx=find(nbasin_faults(l_sr_indx,1)==num_raw(:,81)); 
               mflt_indx=find(nbasin_faults(l_sr_indx,2)==num_fault(:,23));
               flt_area=num_raw(flt_indx,21);
               mflt_area=sum(num_fault(mflt_indx,6).*num_fault(mflt_indx,12));   
            end
            
           elseif id<300
           l_sr_indx=find(num_sec(l_sec_indx,15)==nbasin_faults(:,1));
           else
           l_sr_indx=find(id==nbasin_faults(:,1));    
           end
        end
        
        slip_rate_nb.sr{nbindx}(i)=sslip_rate(i);%store simu slip rate
        
        %replace simulated slip rate with data from lake
         if id>600
            if length(find(id==nbasin_faults(:,2)))>1
            sslip_rate(i)=normrnd(nbasin_faults(l_sr_indx,7),nbasin_faults(l_sr_indx,8));
        %for multifault that includes fault that both do and do not offset reflector
        %include weighted average slip rate of both types of faults
            else
            sslip_rate(i)=((flt_area*normrnd(nbasin_faults(l_sr_indx,3),nbasin_faults(l_sr_indx,4)))+...
              sslip_rate(i)*(mflt_area-flt_area))/mflt_area;
    
            end
         else
            sslip_rate(i)=normrnd(nbasin_faults(l_sr_indx,3),nbasin_faults(l_sr_indx,4));   
         end
        %Renforce postive slip rate
                if sslip_rate(i)<0
                sslip_rate(i)=0.0001;
                end
    end

    if area1<len*((len*1000)^(2/3)*c1(2))/1000 %section is cut off by intersecting fault or thickness of crust
    
    %magnitude and displacement calculated through fault area and c2 (eqs. 7 & 11 of Leonard 2010)   
    sed(i) = 10^(0.5*log10(area1*10^6)+log10(c2_pick(i)));
    
    else
        
    %magnitude and displacement calculated through fault length, c1, and c2 (eqs. 8 & 12 of Leonard 2010) 
    sed(i) = 10^(5/6*log10(len*10^3)+0.5*log10(c1_pick(i))+log10(c2_pick(i)));

    if ~isreal(sed(i))
        sed(i)=0.001;
    end
    
    end
   
    %RI from dis/slip_rate for each simu (e.g. Wallace 1970)
    ri(i)=sed(i)*1000/sslip_rate(i);

end %end simulations for each source

mean_sr(ss) = mean(sslip_rate); std_sr(ss)= std(sslip_rate);
range_sr(ss,:)=[min(sslip_rate) max(sslip_rate)];

%mean ri given sed is log normal dist. and sslip rate is normal dist.
mean_r(ss)=exp(mean(log(sed)))*1000/mean(sslip_rate);

std_r(ss)=fitdist(ri,'Lognormal').sigma;

range_r(ss,:) = [min(ri) max(ri)];

    %Bird2007 test that well constrained slip rate is one
    %with median value > 95% CI
    if median(sslip_rate)>(std_sr(ss)*4)
        constraint_test(ss)=1;
    else
        constraint_test(ss)=0;%poorly constrained
    end
    
end %end simulation for all sources



%% Create matlab variable to store slip rate and reucrrence interval data

MSSD_source_calcs=[MSSD_id(:,1),mean_sr,std_sr,log(mean_r),std_r];
%save('MSSD_source_calcs.mat','MSSD_source_calcs','range_r','range_sr','slip_rate_nb','constraint_test');

% Write data into MSSD spreadsheat

[num_sec,txt_sec,raw_sec]=xlsread('MSSD.xlsx',2);
[num_fault,txt_fault,raw_fault]=xlsread('MSSD.xlsx',3);
[num_mfault,txt_mfault,raw_mfault]=xlsread('MSSD.xlsx',4);

% Make sure spreadsheet is closed when running
% Writes slip rate +/- 1 sigma error
% Writes recurrence interval values +/- 1 sigma error (as error is asymmetric) 

indx=zeros(length(MSSD_source_calcs),1);
simu_count=1;

for i=1:length(MSSD_source_calcs)
    
    if i== 20*simu_count
        disp(['MSSD source: ',num2str(i),' out of ',num2str(length(MSSD_source_calcs))]);
        simu_count = simu_count + 1;
    end
    
    %Write data for fault sections
    if MSSD_source_calcs(i,1)<300
       indx(i)=find(MSSD_source_calcs(i,1)==num_sec(:,1));  
       writematrix(MSSD_source_calcs(i,2:3),'MSSD.xlsx','Sheet','SectionGeometry','Range',['AA',num2str(indx(i)+1)]);
       writematrix([exp(MSSD_source_calcs(i,4)-MSSD_source_calcs(i,5)),exp(MSSD_source_calcs(i,4)), exp(MSSD_source_calcs(i,4)+MSSD_source_calcs(i,5))],...
       'MSSD.xlsx','Sheet','SectionGeometry','Range',['W',num2str(indx(i)+1)]);
    
    elseif MSSD_source_calcs(i,1)>600
       indx(i)=find(MSSD_source_calcs(i,1)==num_mfault(:,1));
       writematrix(MSSD_source_calcs(i,2:3),'MSSD.xlsx','Sheet','MultiFaultGeometry','Range',['F',num2str(indx(i)+1)]);
       writematrix([exp(MSSD_source_calcs(i,4)-MSSD_source_calcs(i,5)),exp(MSSD_source_calcs(i,4)), exp(MSSD_source_calcs(i,4)+MSSD_source_calcs(i,5))],...
       'MSSD.xlsx','Sheet','MultiFaultGeometry','Range',['L',num2str(indx(i)+1)]);
  
    else
      indx(i)=find(MSSD_source_calcs(i,1)==num_fault(:,1)); 
      writematrix(MSSD_source_calcs(i,2:3),'MSSD.xlsx','Sheet','FaultGeometry','Range',['O',num2str(indx(i)+1)]);
      writematrix([exp(MSSD_source_calcs(i,4)-MSSD_source_calcs(i,5)),exp(MSSD_source_calcs(i,4)), exp(MSSD_source_calcs(i,4)+MSSD_source_calcs(i,5))],...
      'MSSD.xlsx','Sheet','FaultGeometry','Range',['T',num2str(indx(i)+1)]);
    end    
end

