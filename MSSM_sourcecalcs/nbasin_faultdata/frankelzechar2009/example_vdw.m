% San Andreas Fault Biskra Palms data
close all;

% AGES REPORTED  BY VAN DER WOERD ET AL 2006 JGR
% ages = [34.6; 33; 36.9; 36.3; 36.8; 36.1; 35.6; 34.2; 39.7; 38.6; 45.8; 35.4; 31; 40.1; 40.5; 33.4; 32.8; 36.4; 37.1; 33.9] * 1000;
% ageUncertainties = [1.6; 4; 1.7; 4.3; 4.4; 5.9; 1.6; 4.4; 4.1; 4.8; 4.2; 3.9; 3.7; 4.3; 2.9; 2.6; 2.6; 2.8; 2.8; 2.7] * 1000;

% "GOOD" AGES REPORTED  BY VAN DER WOERD ET AL 2006 JGR
ages = [34.6; 33; 36.9; 36.3; 36.8; 36.1; 35.6; 34.2; 39.7; 38.6; 35.4; 31; 40.1; 40.5; 33.4; 32.8; 36.4; 37.1; 33.9] * 1000;
ageUncertainties = [1.6; 4; 1.7; 4.3; 4.4; 5.9; 1.6; 4.4; 4.1; 4.8; 3.9; 3.7; 4.3; 2.9; 2.6; 2.6; 2.8; 2.8; 2.7] * 1000;

offsets = [565];
offsetUncertainties = [80];



[t p_t] = ages_gaussian(ages, ageUncertainties, 0.6827);
[d p_d] = displacements_gaussian(offsets, offsetUncertainties, 0.6827);
[v p_v] = slip_rate(t, p_t, d, p_d, 0.6827);

[t p_t] = ages_gaussian(ages, ageUncertainties, 0.9545);
[d p_d] = displacement_gaussian(offsets, offsetUncertainties, 0.9545);
[v p_v] = slip_rate(t, p_t, d, p_d, 0.9545);