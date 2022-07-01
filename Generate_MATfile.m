%CREATE MAT FILES FROM MALAWI AFD SPREADSHEET
% Import data from spreadsheet

clear
close all

[num_sec,txt_sec,raw_sec] = xlsread('MSSM_sourcecalcs/MSSM.xlsx',2);

[num_fault,txt_fault,raw_fault] = xlsread('MSSM_sourcecalcs/MSSM.xlsx',3);

[num_multi_fault,txt_multi_fault,raw_multi_fault] = xlsread('MSSM_sourcecalcs/MSSM.xlsx',4);

[num_MSSM,txt_MSSM,raw_MSSM] = xlsread('MSSM_sourcecalcs/MSSM.xlsx',5);
%Create mat file
save('MSSM_sources')
