%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%   Gutenberg-Richter Relationships for East Africa   %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clear
close all

mydir  = pwd; idcs   = strfind(mydir,'/');
addpath(mydir(1:idcs(end)-1));
load('syncat_PSHA_mssm_input','Region');
addpath([mydir(1:idcs(end)-1) '/gis_files']); addpath([mydir(1:idcs(end)-1) '/misc_functions']);

load EastAfrica_GSHHSdata 
load ISC_GEM_EAcatalog %note this is the same as the GEM catalog


%% Pick and sort data from assessed region around Malawi

Mmin=3;
pick1 = find(EAcatalog(:,7) >= Region(3) & EAcatalog(:,7) <= Region(4) & EAcatalog(:,8) >= Region(1) & EAcatalog(:,8) <= Region(2) & EAcatalog(:,10) >= Mmin);
EAcatalog_R1 = EAcatalog(pick1,:);

%plot map of events
figure 
plot(EAshoreline(:,1),EAshoreline(:,2),'k-'); hold on; 
plot(EAcountry(:,1),EAcountry(:,2),'r-'); hold on;
% plot(EAriver(:,1),EAriver(:,2),'b-'); hold on;
axis equal; axis([10 60 -30 20]);
% plot(EAcatalog_R1(:,8),EAcatalog_R1(:,7),'g.'); hold on;
scatter(EAcatalog_R1(:,8),EAcatalog_R1(:,7),25,EAcatalog_R1(:,10),'filled'); hold on;
plot(Region([1 2 2 1 1]),Region([3 3 4 4 3]),'k-'); hold on;

%% Yearly rate
year(:,1) = 1901:5:2011;
year(:,2) = 1905:5:2015;
year(:,3) = 1902.5:5:2012.5;

for ii = 1:length(year(:,1))
    RateYear_R1(ii,1) = length(find(EAcatalog_R1(:,1) >= year(ii,1) & EAcatalog_R1(:,1) <= year(ii,2)))/5;
end

figure 
plot(year(:,3),RateYear_R1,'b-');
xlabel('Year'); ylabel('Annual event counts');
%%
% Magnitude recurrence rate
YearStart = 1965;

EAcatalog_R1 = EAcatalog_R1(find(EAcatalog_R1(:,1) >= YearStart),:);

mag_check = Mmin:0.1:6.5;
for ii = 1:length(mag_check)
    RateMag_R1(ii,1) = length(find(EAcatalog_R1(:,10) >= mag_check(ii)))/45;
end

% Background seismicity (Zone 9 - Rukwa-Malawi Rift; Poggi et al. 2017) 
b_background    = 1.02;
mmax_background = 7.9;

a_R1 = log10(length(find(EAcatalog_R1(:,10) >= 5))/45) + b_background*4.5;
    
figure
semilogy(mag_check,RateMag_R1(:,1),'rx',mag_check,10.^(a_R1-b_background*mag_check),'b-'); hold on;
semilogy(mmax_background*ones(1,2),[10^-6 10^0],'c--');
xlabel('Magnitude'); ylabel('Annual frequency of exceedance');
axis square; grid on; axis([3 8 10^-4 10^1]);

%Option to save data filted for Malawi into seperate file
%save('ISC_GEM_EAcatalog_Malawi','RateMag_R1','mag_check')

%% Plot Earthquake catalog map for all East Africa

figure
plot(EAshoreline(:,1),EAshoreline(:,2),'k-'); hold on; 
plot(EAcountry(:,1),EAcountry(:,2),'r-'); hold on;
% plot(EAriver(:,1),EAriver(:,2),'b-'); hold on;
axis equal; axis([10 50 -30 20]);
% plot(EAcatalog(:,8),EAcatalog(:,7),'g.'); hold on;

tmp1 = find(EAcatalog(:,10) >= 4 & EAcatalog(:,10) < 5);
tmp2 = find(EAcatalog(:,10) >= 5 & EAcatalog(:,10) < 6);
tmp3 = find(EAcatalog(:,10) >= 6 & EAcatalog(:,10) < 7);
tmp4 = find(EAcatalog(:,10) >= 7);
plot(EAcatalog(tmp1,8),EAcatalog(tmp1,7),'bo','MarkerSize',2); hold on;
plot(EAcatalog(tmp2,8),EAcatalog(tmp2,7),'go','MarkerSize',4); hold on;
plot(EAcatalog(tmp3,8),EAcatalog(tmp3,7),'yo','MarkerSize',6); hold on;
plot(EAcatalog(tmp4,8),EAcatalog(tmp4,7),'ro','MarkerSize',10); hold on;
colormap('jet');

%% Derive 95% CI around instrumental G-R record using Tinti Mulargia 1987 method

Mmax = 6.5; Mc=4.5; %see Hodge et al (2015)

Alpha  = 0.05; K_Alpha = -norminv(Alpha/2,0,1);
mag_check = Mmin:0.2:Mmax; TotalYear=2015-YearStart;

for ii = 1:length(mag_check)
        CumNum(ii) = length(find(EAcatalog_R1(:,10) >= mag_check(ii)));
        AnnualCumNum(ii) = CumNum(ii)/TotalYear;
end

dataindex = find(EAcatalog_R1(:,10) >= Mc);
num_aftershock = length(dataindex);

GRpara(2) = num_aftershock*log10(exp(1))/(sum(EAcatalog_R1(dataindex,10)-Mc+0.05));
GRpara(1) = log10(length(find(EAcatalog_R1(:,10) >= Mc))/TotalYear) + GRpara(2)*Mc;
GRpara(3) = GRpara(2)/sqrt(num_aftershock); % Standard deviation of the estimated b value

% Tinti-Mulargia formula for GR fitting
% MeanMandZ(1:2) = [mean(Data(dataindex,6)) (mean(Data(dataindex,6))-Mc(ijk))/0.1-0.5];
MeanMandZ(1:2) = [mean(EAcatalog_R1(dataindex,10)) (mean(EAcatalog_R1(dataindex,10))-Mc+0.05)/0.1-0.5];
GRparaTM(2) = (-log(MeanMandZ(2)/(MeanMandZ(2)+1))/0.1)/log(10);
GRparaTM(1) = log10(num_aftershock/TotalYear) + GRparaTM(2)*Mc;
GRparaTM(4) = (-log(betaincinv(Alpha/2,MeanMandZ(2)*num_aftershock,num_aftershock))/0.1)/log(10);
GRparaTM(5) = (-log(1- betaincinv(Alpha/2,num_aftershock,MeanMandZ(2)*num_aftershock+1))/0.1)/log(10);    
GRparaTM(3) = (GRparaTM(4)-GRparaTM(5))/(2*K_Alpha);  

semilogy([Mc Mmax],10.^(GRpara(1)-GRpara(2)*[Mc Mmax]),'r-',[Mc Mmax],10.^(GRparaTM(1)-GRparaTM(2)*[Mc Mmax]),'g-','LineWidth',2);
%     title(['ML: a=',num2str(GRpara(1)),' & b=',num2str(GRpara(2)),' vs TM: a=',num2str(GRparaTM(1)),' & b=',num2str(GRparaTM(2))]);
    title(['ML: b=',num2str(GRpara(2)),' & sigmab=',num2str(GRpara(3)),' vs TM: b=',num2str(GRparaTM(2)),' & sigmab=',num2str(GRparaTM(3))]);
    legend(num2str(num_aftershock));
    
    
%tmp1 = find(EAcatalog_R1(:,10) >= 4 & EAcatalog_R1(:,10) < 5);
%tmp2 = find(EAcatalog_R1(:,10) >= 5 & EAcatalog_R1(:,10) < 6);
%tmp3 = find(EAcatalog_R1(:,10) >= 6 & EAcatalog_R1(:,10) < 7);

% Confidence interval of the GR relationship for the entire region
CompleteYear = [5.0 5.5 YearStart; 5.5 6.0 YearStart; 6.0 6.5 YearStart];

for kk = 1:length(CompleteYear)
    tmp2 = find(EAcatalog_R1(:,10) >= CompleteYear(kk,1) & EAcatalog_R1(:,10) < CompleteYear(kk,2) & CompleteYear(kk,3) <= EAcatalog_R1(:,1)); 
    CountW(kk) = length(tmp2);
    TimeW(kk)  = 2015 - CompleteYear(kk,3) + 1;
    %clear tmp2

end

GRFitW = GRrelation_MLEWeichert_EQMATca(CountW,5.25:0.5:6.25,0.5,TimeW);

num_sample = 10000;
N0 = exp((log(GRFitW(3)) - 0.5*log(1+(GRFitW(4)/GRFitW(3))^2)) + sqrt(log(1+(GRFitW(4)/GRFitW(3))^2))*normrnd(0,1,num_sample,1));
b  = GRFitW(1) + GRFitW(2)*normrnd(0,1,num_sample,1);

log10N = log10(N0)*ones(1,42) - b*(4:0.1:8.1);

for kk = 1:42;
      Stat_log10N(kk,1:7) = prctile(log10N(:,kk),[5 10 16 50 84 90 95]);
end

figure (30);
semilogy(mag_check,AnnualCumNum,'bo',5.25:0.5:6.25,CountW./TimeW,'rs',4:0.1:8.1,10.^(Stat_log10N(:,4)),'c-',4:0.1:8.1,10.^(Stat_log10N(:,[3 5])),'c--','LineWidth',1);
        axis square; xlabel('Magnitude'); ylabel('N(>m)'); axis([Mmin Mmax 0.01 100]); hold on;

