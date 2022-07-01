close all;

% % % % % % % % % % 
% S. Tenjindo data
% % % % % % % % % %

% interpret displacement range as boxcar
[t p_t] = age_boxcar(60000, 70000, 0.6827, false, [], []);
[d p_d] = displacements_boxcar(70, 90, 0.6827, false, [], []);
[v p_v] = slip_rate(t, p_t, d, p_d, 0.6827, true, 0.6, 2);

[t p_t] = age_boxcar(60000, 70000, 0.9545, false, [], []);
[d p_d] = displacements_boxcar(70, 90, 0.9545, false, [], []);
[v p_v] = slip_rate(t, p_t, d, p_d, 0.9545, false, [], []);

% interpret displacement range as Gaussian 1 std
[t p_t] = age_boxcar(60000, 70000, 0.6827, false, [], []);
[d p_d] = displacements_gaussian(80, 10, 0.6827, false, [], []);
[v p_v] = slip_rate(t, p_t, d, p_d, 0.6827, true, 0.6, 2);

[t p_t] = age_boxcar(60000, 70000, 0.9545, false, [], []);
[d p_d] = displacements_gaussian(80, 10, 0.9545, false, [], []);
[v p_v] = slip_rate(t, p_t, d, p_d, 0.9545, false, [], []);