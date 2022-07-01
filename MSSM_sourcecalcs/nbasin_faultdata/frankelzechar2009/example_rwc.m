% Red Wall Canyon data
close all;

% DISPLACEMENT IS TREATED AS SUM OF GAUSSIANS, ONE FOR EACH STREAM OFFSET; 
% MEAN TAKEN TO BE THALWEG OFFSET, STANDARD DEVIATION TAKEN TO BE
% 1/2 THE CHANNEL WIDTH
disps = [307; 291; 304; 286; 301; 306; 287]; 
dispUncertainties = [24; 13; 12; 21; 17; 22; 33] / 2;

% "GOOD" AGES USED IN KURT'S PAPER
ages = [76900; 62500; 70000; 65200; 70800; 86400; 69900; 69400; 72100; 55700; 81900; 64800; 62300]; 
ageUncertainties = [1800; 2400; 1800; 2700; 2500; 2400; 2100; 2100; 2300; 2000; 2500; 2500; 2100]; 

[t p_t] = ages_gaussian(ages, ageUncertainties, 0.6827, true, [], []);
[d p_d] = displacements_gaussian(disps, dispUncertainties, 0.6827, true, [], []);
[v p_v] = slip_rate(t, p_t, d, p_d, 0.6827, true, [], []);
%%
[t p_t] = ages_gaussian(ages, ageUncertainties, 0.9545, false, 0, 0);
[d p_d] = displacements_gaussian(disps, dispUncertainties, 0.9545, false, 0, 0);
[v p_v] = slip_rate(t, p_t, d, p_d, 0.9545, false, 0, 0);
%%
% ALL AGES REPORTED IN KURT'S PAPER
ages = [219400; 76900; 124000; 62500; 70000; 65200; 70800; 86400; 69900; 69400; 72100; 55700; 81900; 64800; 62300; 36500]; 
ageUncertainties = [3500; 1800; 2400; 2400; 1800; 2700; 2500; 2400; 2100; 2100; 2300; 2000; 2500; 2500; 2100; 5600]; 

[t p_t] = ages_gaussian(ages, ageUncertainties, 0.6827, true, [], []);
[d p_d] = displacements_gaussian(disps, dispUncertainties, 0.6827, false, 0, 0);
[v p_v] = slip_rate(t, p_t, d, p_d, 0.6827, true, [], []);

[t p_t] = ages_gaussian(ages, ageUncertainties, 0.9545, false, 0, 0);
[d p_d] = displacements_gaussian(disps, dispUncertainties, 0.9545, false, 0, 0);
[v p_v] = slip_rate(t, p_t, d, p_d, 0.9545, false, 0, 0);