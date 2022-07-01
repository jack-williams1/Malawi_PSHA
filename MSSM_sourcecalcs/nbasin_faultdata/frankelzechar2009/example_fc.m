% Furnace Creek data
close all;

% ALL AGES REPORTED BY FRANKEL 2007 (INCL. OUTLIER)
ages = [100.7; 111.6; 64.3; 78.7; 83.8; 93.5; 86.4; 103.7; 97.4] * 1000; 
ageUncertainties = [8.3; 8.4; 4.7; 7.8; 8.4; 7.1; 6.5; 7.6; 7.3] * 1000; 

[t p_t] = ages_gaussian(ages, ageUncertainties, 0.6827, true, 20000, 160000);
[d p_d] = displacements_gaussian(290, 20, 0.6827, false, 0, 0);
[v p_v] = slip_rate(t, p_t, d, p_d, 0.6827, true, 1, 7);

[t p_t] = ages_gaussian(ages, ageUncertainties, 0.9545, false, 0, 0);
[d p_d] = displacements_gaussian(290, 20, 0.9545, false, 0, 0);
[v p_v] = slip_rate(t, p_t, d, p_d, 0.9545, false, 0, 0);


% "GOOD" AGES USED BY FRANKEL ET AL 2007 GRL
ages = [100.7; 111.6; 78.7; 83.8; 93.5; 86.4; 103.7; 97.4] * 1000; 
ageUncertainties = [8.3; 8.4; 7.8; 8.4; 7.1; 6.5; 7.6; 7.3] * 1000; 

[t p_t] = ages_gaussian(ages, ageUncertainties, 0.6827, true, 20000, 160000);
[d p_d] = displacements_gaussian(290, 20, 0.6827, false, 0, 0);
[v p_v] = slip_rate(t, p_t, d, p_d, 0.6827, true, 1, 7);

[t p_t] = ages_gaussian(ages, ageUncertainties, 0.9545, false, 0, 0);
[d p_d] = displacements_gaussian(290, 20, 0.9545, false, 0, 0);
[v p_v] = slip_rate(t, p_t, d, p_d, 0.9545, false, 0, 0);