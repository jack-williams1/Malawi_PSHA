% DETERMINES HANGING WALL FLEXURE IN NORMAL FAULTS FROM EQUATION S1 OF MUIRHEAD ET AL 2016
% BASED ON EQUATIONS DETAILED IN BILLINGS AND KATTENHRON (2005)

% INPUTS FOR ALPHA ARE:
% Te: Thickness of elastic crust (m)
% E: Young Modulus
% v: Poissons Ratio
% p: Density (m/kg^3)
% g: Gravity (m/s^2)

% Wo: Max hanging-wall deflection per basin (m)
% b_width: Basin width (km)

%% Input and Sort Data
num_table=readtable('HangingWallFlexureInputs.xlsx','Sheet','inputs','Range','A3:O11');%Update cell range if add more basins. Currently set for south Malawi only
%[txt,raw]=xlsread('HangingWallFlexureInputs','inputs','A3:A11');% Include basin names for plots

num_data=table2array(num_table(:,2:width(num_table)));

Te=[num_data(:,1)-num_data(:,2) num_data(:,1) num_data(:,1)+num_data(:,2)];
E=[num_data(:,3)-num_data(:,4) num_data(:,3) num_data(:,3)+num_data(:,4)];
v=num_data(:,5);
p=[num_data(:,6)+num_data(:,7) num_data(:,6) num_data(:,6)-num_data(:,7)];
g=num_data(:,8);

Wo = [num_data(:,12)-num_data(:,13) num_data(:,12) num_data(:,12)+num_data(:,13)];
b_width = num_data(:,14);

%% Calculate Alpha Parameter from Eq S2 of Muirhead et al (2016)
% Varies for each basin and elastic properties of crust

a=zeros(length(Te),3);
for i =1:length(Te)
    for j=1:3
    a(i,j)=((E(i,j)*Te(i,j)^3)/(3*p(i,j)*g(i)*(1-v(i)^2)))^(0.25);
    end
end

%% Calculate and hanging wall flexure profile (W) over 100 km distance (in metres)

distance=[0:10^3:10^5]; %Create sequence of intergers from 0-100

%Maximum deflection requires highest Wo (i.e. T_throw) and lowest a (i.e. Te)
a=fliplr(a);

for i=1:length(Te)%Number of basins
      for j=1:3 
            for k=1:length(distance)
                deflection(i,j,k)=Wo(i,j)*(exp(-distance(k)/a(i,j)))*(cos(distance(k)/a(i,j))); %Calculates magnitude of deflection from Eq S1
                profile(i,j,k)=(Wo(i,j)-deflection(i,j,k))-Wo(i,j); %Calculates hanging wall profile after deflection
        end
    end
end

%% Determine strain profile using equation S4 of Muirhead et al (2016)
syms xx % Create symbol xx
strainCalc=cell(3,length(Te));

for i=1:length(Te)
    for j=1:3
        deflectionCalc(i,j)=Wo(i,j)*exp(-xx/a(i,j))*(cos(xx/a(i,j))); %Equation S1 of Muirhead et al (2016) in terms of xx
        derv(i,j)=diff(diff(deflectionCalc(i,j)),xx); %Determines second deravative of equation S1
        strainCalc{i,j}=derv(i,j);%Place each differentation equation into a cell array
        strainCalc_{i,j}=symfun(strainCalc{i,j},xx);%Turn equations from cell array into functions
    end
end


for i=1:length(Te)
    for j=1:3
        for k=1:length(distance)
            strain(i,j,k)=strainCalc_{i,j}(distance(k))*Te(i,j)/2;% Run equation S4 of Muirhead et al (2016) where Y is half thickness of plate
            strain_pc(i,j,k)=(double(strain(i,j,k)))*100; % Convert strain data type to double and %
  
        end
    %Average extensional strain over region that spans basin width
          meanExtensionStrain(i,j)=mean(strain_pc(i,j,(1:b_width(i))));
    %Total extensional displacement acriss basin
           totalExtensionDisplacement(i,j)=b_width(i)*1000*meanExtensionStrain(i,j)/100;% %Ignore for Lower Shire where different widths
   end
end

%Note meanExtension and totalExtension for Makanjira grabens calculated
%below, once contribution from each border fault is considered

%% Sort data for plots

%For plotting, convert distance to km
for i=1:length(distance)
    distance_K(i)=distance(i)/1000;
end

% CONVERT MAKANJIRA RESULTS SO THEY CAN BE SHOWN ON SAME PLOT ACROSS 90 KM WIDE REGION

mak_e_idx=find(string(table2cell(num_table(:,1)))=='Makanjira-East');
mak_w_idx=find(string(table2cell(num_table(:,1)))=='Makanjira-West');
%%

mak_totProfile=zeros(3,b_width(mak_e_idx));
mak_totStrain=zeros(3,b_width(mak_e_idx));
w_sp=100-b_width(mak_e_idx)+1;

for i=1:3
 %Flip arrays for Makanjira East so runs west to east
mak_e_Profile(i,:)=flip(profile(mak_e_idx,i,:));
mak_w_Profile(i,:)=profile(mak_w_idx,i,:);
mak_e_strain_pc(i,:)=flip(strain_pc(mak_e_idx,i,:));
mak_w_strain_pc(i,:)=strain_pc(mak_w_idx,i,:);
%Combine with Makanjira Profile East across 90 km rift

mak_totProfile(i,:) =  mak_e_Profile(i,1:b_width(mak_e_idx))+mak_w_Profile(i,w_sp:100);
mak_totStrain_pc(i,:) =  mak_e_strain_pc(i,1:b_width(mak_e_idx))+mak_w_strain_pc(i,w_sp:100);

mak_meanExtensionStrain(i,:)=mean(mak_totStrain_pc(i,:));%Revised strain for Makanjira
mak_totalExtension(i,:)=(90*1000*mak_meanExtensionStrain(i,:)/100); %Revised extension for Makanjira
end
% Add new values to mean Extension Strain matrix
meanExtensionStrain(mak_e_idx,:)=mak_meanExtensionStrain';
meanExtensionStrain(mak_w_idx,:)=mak_meanExtensionStrain';
totalExtensionDisplacement(mak_e_idx,:)=mak_totalExtension';
totalExtensionDisplacement(mak_w_idx,:)=mak_totalExtension';

%% Calculate strain at points along profile where fault

num_table1=readtable('HangingWallFlexureInputs.xlsx','Sheet','dist2bf','Range','B2:D70');
faultstrain=zeros(length(table2cell(num_table1)),3);

for i=1:length(faultstrain)
    b_indx=find(string(table2cell(num_table1(i,1)))==table2cell(num_table(:,1)));%index fault basin in num_table
    
    if table2array(num_table1(i,3))>0 %Correction for Cassimo fault which is behind border fault
        b_dist_indx=find(round(table2array(num_table1(i,3)))==distance_K);%index fault distance along profile

    if b_indx==mak_e_idx
        faultstrain(i,:)=mak_totStrain_pc(:,b_dist_indx);
    else
        faultstrain(i,:)=strain_pc(b_indx,:,b_dist_indx);
    end
    
    end
end

%% Save variables for further analysis and plotting
save('hangingwallflexcalcs','strain_pc','profile','meanExtensionStrain','totalExtensionDisplacement','deflection',...
    'mak_totProfile','mak_totStrain_pc','num_table','distance','b_width');

% save mean strain values in spreadsheet
writematrix(meanExtensionStrain,'HangingWallFlexureInputs.xlsx','Sheet','inputs','Range','W3:Y12');
writematrix(totalExtensionDisplacement,'HangingWallFlexureInputs.xlsx','Sheet','inputs','Range','Z3:AB12');

writematrix(faultstrain,'HangingWallFlexureInputs.xlsx','Sheet','dist2bf','Range','E2:G70');


