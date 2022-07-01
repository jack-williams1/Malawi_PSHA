%% Derive slip rates for NB Faults 

% Use data from Donna Shillington sent 29/03/21
% Model uncertainity from age of reflector and fault dip
% using methods of frankel and zechar 2009

%Input and sort data
addpath('NorthBasinFaultSlipRates')
addpath('frankelzechar2009')

[num_xref,txt_xref,raw_xref]=xlsread('nbasin_faults_idxref.xlsx');

%selection of sd error (1 or 2)
sd_sel = 1;

if sd_sel ==1
sd_opt = 0.6827;
elseif sd_sel==2
sd_opt = 0.9545;
end

sr_est=zeros(16,3);

%range of fault dips considered
f_dip = [65 53 40];

%% Slip rate calculations

%Uncertainity in reflector age (Scholz et al 2007)
%Type of error not report, but assume best represented by one sigma
[reflector_t reflector_p_t] = ages_gaussian(74200, 5300, sd_opt, false, [], []);

v_err = 7.5;%vertical error of throw measurement;

for i=1:16

    %Read in table
    %Note file name for fault 16 amended to allow easy reading of data
    %Also delimiters revised to tabs
    faultdata{i}=readtable(['fault',num2str(i),'_offsets_feb2019_megadrought_ll.txt'],...
        'ReadVariableNames', false,'Delimiter','\t');
    
    fault_name(i,1)=string(['fault ',num2str(i)]);
    
    %extract all non-zero offsets
    if all(table2array(faultdata{i}(:,3)))==1
        offset{i} = table2array(faultdata{i}(:,3));
    else
        tmp_indx=find(table2array(faultdata{i}(:,3))~=0);
        offset{i} = table2array(faultdata{i}(tmp_indx,3)); clear tmp_indx
    end
    
    %min displacement estimate from steepest fault dip
    disp_min = (offset{i}-v_err)/sin(f_dip(1)*pi/180); 
    %max displacement estimate from lowest fault dip
    disp_max = (offset{i}+v_err)/sin(f_dip(3)*pi/180);
    
    %obtain prob density function for offset measurements
    %where >1 offset measurement, prob density taken for each measurement
    %and then a cumulative prob density derived from a normalisation
    
    [d_fault{i} pd_fault{i}]=displacements_boxcar(disp_min,disp_max,sd_opt,false,[],[]);
    
    %derive slip rate for each fault
    [v{i},pd_v{i},sr_est(i,1),sr_est(i,2),sr_est(i,3)]=slip_rate(reflector_t, reflector_p_t, d_fault{i}, pd_fault{i}, sd_opt, false, [], []);
    
    clear disp disp_sd

end

%% Read in data from MSSM to calc multifault slip rates from weighted area averages

mydir  = pwd; idcs   = strfind(mydir,'/');
addpath(mydir(1:idcs(end)-1));

num_fault=readtable('MSSM.xlsx','Sheet','FaultGeometry');
fault_area=readtable('MSSM.xlsx','Sheet','Leonard2010','Range','A:Y');

nb_multifault_xref=table2array(num_fault(:,[1,23]));
area_xref=table2array(fault_area(:,1));
mflt_indx1=zeros(height(num_xref),8);

for i=1:length(mflt_indx1)
    
    if num_xref(i,3)>0

    mflt_indx1(i,5)=find(num_xref(i,3)==nb_multifault_xref(:,1)); %index in 'fault geom' tab
    mflt_indx1(i,6)=find(num_xref(i,3)==area_xref(:,1),1);%index in 'Leonard 2010' tab
    end
 
end

indxclr=find(mflt_indx1(:,5)==0);%clear faults not mapped in the MSSM
mflt_indx1(indxclr,:)=[]; 
mflt_indx1(:,1:2)=[nb_multifault_xref(mflt_indx1(:,5),1),nb_multifault_xref(mflt_indx1(:,5),2)];

%add in slip rates from sr-est
for i=1:height(mflt_indx1)
    tmp=find(mflt_indx1(i,1)==num_xref(:,3));
    mflt_indx1(i,3)=sr_est(num_xref(tmp,1),1);
    mflt_indx1(i,4)=sr_est(num_xref(tmp,1),1)-sr_est(num_xref(tmp,1),2);
    clear tmp
end
%% Calc weighted average slip rates

mflt_id=unique(mflt_indx1(:,2));

for i=1:length(mflt_id)

if mflt_id(i)~=0 %ignore non multifault faults
    
    mflt_indx2=find(mflt_id(i)==mflt_indx1(:,2));
    
    for ii=1:length(mflt_indx2)
       %find area, slip rate and 1 sigma error of each participating fault
       f_area(ii)=fault_area.Var22(mflt_indx1(mflt_indx2(ii),6));
       sr_mean(ii)=f_area(ii)*mflt_indx1(mflt_indx2(ii),3);
       sr_std(ii)=f_area(ii)*mflt_indx1(mflt_indx2(ii),4);
    end
    
   %calc area weighted average of multifault slip rates 
   mflt_indx1(mflt_indx2,7)=sum(sr_mean)/sum(f_area);
   mflt_indx1(mflt_indx2,8)=sum(sr_std)/sum(f_area);
   clear mflt_indx2 f_area sr_mean sr_std
end

end
%% Save outputs

%create table and save as csv file for all NB faults
nb_slipratetable=table(fault_name,sr_est);
writetable(nb_slipratetable,'nb_slipratetable.csv')

%create table for NB faults included in the MSSM only
writetable(table(mflt_indx1),'nb_slipratetable_id.csv');

%save slip rate values for plotting in MSSD source plots
save('v_id.mat','v');