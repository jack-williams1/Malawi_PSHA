% Indian Creek data
close all;
ages = [63.3; 66.8; 74.5; 81.2; 59.2; 77.4; 73.8; 74.9] * 1000;
ageUncertainties = [4.7; 5; 5.6; 6.1; 4.4; 5.8; 5.5; 6] * 1000;

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% Get 4 stream offset measurements and errors and use them here
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
disps = [172; 174; 188; 179]; 
dispUncertainties = [8; 11; 16; 38] / 2;

[t p_t] = ages_gaussian(ages, ageUncertainties, 0.6827, true, 30000, 120000);
[d p_d] = displacements_gaussian(disps, dispUncertainties, 0.6827, true, 80, 280);
[v p_v] = slip_rate(t, p_t, d, p_d, 0.6827, true, 1, 5);

[t p_t] = ages_gaussian(ages, ageUncertainties, 0.9545, false, 0, 0);
[d p_d] = displacements_gaussian(disps, dispUncertainties, 0.9545, false, 0, 0);
[v p_v] = slip_rate(t, p_t, d, p_d, 0.9545, false, 0, 0);