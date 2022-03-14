%CREATE MAT FILES FROM MALAWI AFD SPREADSHEET
% Import data from spreadsheet

clear
close all

[num_sec,txt_sec,raw_sec] = xlsread('MSSD_sourcecalcs/MSSD.xlsx',2);

[num_fault,txt_fault,raw_fault] = xlsread('MSSD_sourcecalcs/MSSD.xlsx',3);

[num_multi_fault,txt_multi_fault,raw_multi_fault] = xlsread('MSSD_sourcecalcs/MSSD.xlsx',4);

[num_MSSD,txt_MSSD,raw_MSSD] = xlsread('MSSD_sourcecalcs/MSSD.xlsx',5);
%Create mat file
save('MSSD_sources')
