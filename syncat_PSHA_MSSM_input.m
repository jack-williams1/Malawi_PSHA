%% INPUT SYNCAT_PSHA_MSSM FROM EXCEL TO MATLAB %%

clear all

%READS IN 'syncat_PSHA_MSSM_input.xlsx' AND CONVERTS ALL INPUTS TO A MATLAB VARIABLE

inputexcel = 'syncat_PSHA_MSSM_input.xlsx';

rnd_gen_index = xlsread(inputexcel,1,'C4');
rng(rnd_gen_index,'twister');

%MSSD Source Calc Inputs

dip_w = [xlsread(inputexcel,1,'C21') xlsread(inputexcel,1,'D21') xlsread(inputexcel,1,'E21')]; %dip weighting  
hwf_w = [xlsread(inputexcel,1,'C23') xlsread(inputexcel,1,'D23') xlsread(inputexcel,1,'E23')]; %flexural correction weigting
ext_w = [xlsread(inputexcel,1,'C22') xlsread(inputexcel,1,'D22') xlsread(inputexcel,1,'E22')]; %extension correction weighting

%Syncat options

%Load in name for MSSD syncat
[~,syncat_name] = xlsread(inputexcel,1,'C31:C31');

%Length of syncat
t_limit = xlsread(inputexcel,1,'C19:C19'); 
NumSimu =xlsread(inputexcel,1,'C18');

bg_mmin = xlsread(inputexcel,1,'C43');

%Rupture weightings (Direct MSSD)
weight_SWM(1,1:3) = xlsread(inputexcel,1,'C8:E8');
weight_SWM(2,1:3) = 1:3;

%Fault width weights (Adapted MSSD)
fault_width_weight = [xlsread(inputexcel,1,'C13') xlsread(inputexcel,1,'D13')];
%Recurrence type weights (Adapted MSSD)
recurrence_type_weight = [xlsread(inputexcel,1,'C14') xlsread(inputexcel,1,'D14')];

%Magntiude uncertainity
para_Mw = xlsread(inputexcel,1,'C9:E9'); 

% Depth uncertainty: multiple discrete values
num_TD    = xlsread(inputexcel,1,'C10:C10');
weight_TD = xlsread(inputexcel,1,'C11:L12');
weight_TD = weight_TD(1:2,1:num_TD);


%GMPE options
% 1)  Boore & Atkinson (2008) for shallow crustal events (NGA-WEST1)
% 2)  Campbell & Bozorgnia (2008) for shallow crustal events (NGA-WEST1)
% 3)  Chiou & Youngs (2008) for shallow crustal events (NGA-WEST1)
% 4)  Akkar & Bommer (2010) - Mw = 5.0-7.6; constant sigma model
% 5)  Boore, Stewart, Seyhan & Atkinson (2014) for shallow crustal events (NGA-WEST2)
% 6)  Akkar, Sandikkaya, & Bommer (2014) for shallow crustal events (Rjb model)
% 7)  Akkar, Sandikkaya, & Bommer (2014) for shallow crustal events (Repi model)
% 8)  Akkar, Sandikkaya, & Bommer (2014) for shallow crustal events (Rhypo model)
% 9)  Cauzzi, Faccioli, Vanini, & Bianchini (2015) for shallow crustal events
% 10) Chiou & Youngs (2014) for shallow crustal events
% 11) Atkinson and Adams (2013) for stable continental events (6th Seismic Hazard Model of Canada) - best/lower/upper branches 
% 12) Goulet et al. (2017) for stable continental events (NGA-East preliminary with site amplification; 6th Seismic Hazard Model of Canada) - models 1 to 13 

GMPEindex        = xlsread(inputexcel,1,'C36:L36');% Choses GMPE from excel input file
num_GMPEindex    = length(GMPEindex);
weight_GMPEindex = xlsread(inputexcel,1,'C37:L37');
FMtype           = xlsread(inputexcel,1,'C38:D38');
FMweight         = xlsread(inputexcel,1,'C39:D39');     
Z10_Z25          = xlsread(inputexcel,1,'C40:D40');%Depths to VS = 1 km/s and VS = 2.5 km/s

vs30_site_ref = xlsread(inputexcel,1,'D41');%Reference value for vs30 also to be evaluated

if xlsread(inputexcel,1,'C41') ~= 0
    vs30_site = xlsread(inputexcel,1,'C41')*ones(num_site,1); 
end

%PSHA options

par_opts = [1:1:16];%number of fake parallelisations in PSHA runs on BlueCrystal

obs_duration        = NumSimu*t_limit;
lower_return_period = 100;
num_need            = obs_duration/lower_return_period;%i.e. 10000

Region =  [xlsread(inputexcel,1,'D47') xlsread(inputexcel,1,'F47') xlsread(inputexcel,1,'C47') xlsread(inputexcel,1,'E47')];% For all Malawi. Has to be increments of grid interval
PSHA_Zone = [xlsread(inputexcel,1,'C47') xlsread(inputexcel,1,'D47'); xlsread(inputexcel,1,'C47') xlsread(inputexcel,1,'F47'); xlsread(inputexcel,1,'E47') xlsread(inputexcel,1,'F47'); xlsread(inputexcel,1,'E47') xlsread(inputexcel,1,'D47'); xlsread(inputexcel,1,'C47') xlsread(inputexcel,1,'D47')];

% Grid locations for seismic hazard maps
site_corner_SW     = xlsread(inputexcel,1,'C47:D47');
site_corner_NE     = xlsread(inputexcel,1,'E47:F47');
site_grid_interval = xlsread(inputexcel,1,'C48'); %Grid interval set in spreadsheet at 0.2 degrees

lower_return_period = xlsread(inputexcel,1,'C27');
figure_option = xlsread(inputexcel,1,'C52');

%Max areal source event size in Malawi
max_bg = xlsread(inputexcel,1,'D43');

%Catalog options
cat_option = xlsread(inputexcel,1,'C33');
cat_w = [xlsread(inputexcel,1,'C32') xlsread(inputexcel,1,'D32')];  

% Site options: Coordinates of major locations in Malawi for site_specific PSHA

site_opt=[-13.984404 33.784203;... % Lilongwe (1)
-13.777999 34.459123;... % Salima (2)
-11.460399 34.015926;... % Mzuzu (3)
-9.933439 33.933163;... % Karonga (4)
-15.385787 35.333348;... % Zomba (5)
-15.787254 35.006156;... % Blantyre (6)
-15.065890 35.233263;... % Liwonde (7)
-14.986249 34.956239;... % Balaka (8)
-14.483827 35.257881;... % Mangochi (9)
-14.429176 34.597520;... % Golomoti (10)
-14.221951 34.514127;... % Mtakataka (11)
-16.916607 35.266689]; % Nsanje (12)

site_name_opt=strcat(["Lilongwe","Salima","Mzuzu","Karonga","Zomba","Blantyre","Liwonde",...
    "Balaka","Mangochi","Golomoti","Mtakataka","Nsanje"]); 

site=site_opt(xlsread(inputexcel,1,'C29:E29'),:);
site_name=site_name_opt(xlsread(inputexcel,1,'C29:E29'));

% Spectral vibration periods; "T = 0" indicates the PGA.
T_map = xlsread(inputexcel,1,'C50');% Spectral vibration periods for PSHA Map; 
T =  xlsread(inputexcel,1,'C35:J35'); %for site specific PSHA
num_T = length(T);

group_size=xlsread(inputexcel,1,'C28');
if group_size<num_need
    print('check')
end

prob_level = [xlsread(inputexcel,1,'C42') xlsread(inputexcel,1,'D42')];

%% Save variables under Matlab variable 

save('syncat_PSHA_MSSM_input');