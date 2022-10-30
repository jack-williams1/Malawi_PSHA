%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%   Testing for effect of filtering out Karonga  %%%%%
%%%%%    earthquakes on areal source seismicity  %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
close all

mydir  = pwd; idcs   = strfind(mydir,'/');
addpath(mydir(1:idcs(end)-1));
load('syncat_PSHA_mssm_input','Region','NumSimu','t_limit');
addpath([mydir(1:idcs(end)-1) '/gis_files']); addpath([mydir(1:idcs(end)-1) '/misc_functions']);

load EastAfrica_GSHHSdata 
load ISC_GEM_EAcatalog %note this is the same as the GEM catalog
load new_SZ %with updated source zones for Mozambique
load syncat_bg

%% Pick and sort data from Rukwa-Malawi Source Zone

Mmin=3;
SZ_poly9=polyshape(new_SZ(9).ZONE(:,1),(new_SZ(9).ZONE(:,2))); % Malawi Rift

flagIsIn = inpolygon(EAcatalog(:,8),EAcatalog(:,7),SZ_poly9.Vertices(:,1),SZ_poly9.Vertices(:,2));  
pick1=find(flagIsIn==1);%Filter 1, select only events in Source Zone 9
EAcatalog_R1 = EAcatalog(pick1,:);

pick2 = find(EAcatalog_R1(:,10) >= Mmin);%Filter 2, remove events M<3
EAcatalog_R1 = EAcatalog_R1(pick2,:);

% karonga earthquakes index (RUN ONLY ONCE)
KarongaIndex=[166, 183, 186, 190]; 

EAcatalog_R1_filtered=EAcatalog_R1;
EAcatalog_R1_filtered(KarongaIndex,:)=[];

%% Derive 95% CI around instrumental G-R record using Tinti Mulargia 1987 method
% Plots G-R relationship for no filtering and filtering out Karonga earthquakes 

Mmax = 7; Mc=4.5; YearStart=1965; %see Hodge et al (2015)

Alpha  = 0.05; K_Alpha = -norminv(Alpha/2,0,1);
mag_check = Mmin:0.2:Mmax; TotalYear=2015-YearStart;

%sort syncat_bg for events in SZ 9 and M<7
bg_index=find(bg_catalog(:,6)==3 & bg_catalog(:,4)<=Mmax);
bg_catalog=bg_catalog(bg_index,:);

CumNum=zeros(length(mag_check),1); AnnualCumNum=zeros(length(mag_check),1);


for ii = 1:length(mag_check)
        CumNum(ii) = length(find(EAcatalog_R1(:,10) >= mag_check(ii) & EAcatalog_R1(:,1) >= YearStart));
        AnnualCumNum(ii) = CumNum(ii)/TotalYear;
end

%allow a finer discretization for the bg_syncat
mag_check_syncat = (Mc:0.05:Mmax); bg_AnnualRate=zeros(length(mag_check_syncat),1);

for ii=1:length(mag_check_syncat)
    
    bg_AnnualRate(ii) = length(find(bg_catalog(:,4) >= mag_check_syncat(ii)))/(NumSimu*t_limit); %Rate in Malawi from bg catalog

end

%Plots MFD for SSA-GEM catalog for the SSA GEM catalog
%and for when the catalog is filtered for the 2009 Karonga events

figure (30);

for ii=1:2

    if ii==1 
        tmp_catalog=EAcatalog_R1;%plot MFD for SSA-GEM catalog with Karonga earthquakes
        
    else   
        tmp_catalog=EAcatalog_R1_filtered;%plot MFD for SSA-GEM catalog without Karonga earthquakes
    
    
        for ff = 1:length(mag_check)
            CumNum(ff) = length(find(tmp_catalog(:,10) >= mag_check(ff) & tmp_catalog(:,1) >= YearStart));
            AnnualCumNum(ff) = CumNum(ff)/TotalYear;
        end
  
    end

    dataindex = find(tmp_catalog(:,10) >= Mc);
    num_aftershock = length(dataindex);
    
    MeanMandZ(1:2) = [mean(tmp_catalog(dataindex,10)) (mean(tmp_catalog(dataindex,10))-Mc+0.05)/0.1-0.5];
    GRpara(2) = num_aftershock*log10(exp(1))/(sum(tmp_catalog(dataindex,10)-Mc+0.05));
    GRpara(1) = log10(length(find(tmp_catalog(:,10) >= Mc))/TotalYear) + GRpara(2)*Mc;
    GRpara(3) = GRpara(2)/sqrt(num_aftershock); % Standard deviation of the estimated b value

    % Tinti-Mulargia formula for GR fitting
    % MeanMandZ(1:2) = [mean(Data(dataindex,6)) (mean(Data(dataindex,6))-Mc(ijk))/0.1-0.5];
    
    GRparaTM(2) = (-log(MeanMandZ(2)/(MeanMandZ(2)+1))/0.1)/log(10);
    GRparaTM(1) = log10(num_aftershock/TotalYear) + GRparaTM(2)*Mc;
    GRparaTM(4) = (-log(betaincinv(Alpha/2,MeanMandZ(2)*num_aftershock,num_aftershock))/0.1)/log(10);
    GRparaTM(5) = (-log(1- betaincinv(Alpha/2,num_aftershock,MeanMandZ(2)*num_aftershock+1))/0.1)/log(10);    
    GRparaTM(3) = (GRparaTM(4)-GRparaTM(5))/(2*K_Alpha);  


    % Confidence interval of the GR relationship for the entire region
    CompleteYear = [5.0 5.5 YearStart; 5.5 6.0 YearStart; 6.0 6.5 YearStart];

    for kk = 1:length(CompleteYear)
        tmp2 = find(tmp_catalog(:,10) >= CompleteYear(kk,1) & tmp_catalog(:,10) < CompleteYear(kk,2) & CompleteYear(kk,3) <= tmp_catalog(:,1)); 
        CountW(kk) = length(tmp2);
        TimeW(kk)  = 2015 - CompleteYear(kk,3) + 1;
        %clear tmp2
    end
    
    %Run Tinti and Mulargia function
    GRFitW = GRrelation_MLEWeichert_EQMATca(CountW,5.25:0.5:6.25,0.5,TimeW);

    num_sample = 10000;
    N0 = exp((log(GRFitW(3)) - 0.5*log(1+(GRFitW(4)/GRFitW(3))^2)) + sqrt(log(1+(GRFitW(4)/GRFitW(3))^2))*normrnd(0,1,num_sample,1));
    b  = GRFitW(1) + GRFitW(2)*normrnd(0,1,num_sample,1);

    log10N = log10(N0)*ones(1,42) - b*(4:0.1:8.1);

    for kk = 1:42
        Stat_log10N(kk,1:7) = prctile(log10N(:,kk),[5 10 16 50 84 90 95]);
    end
    
    
    if ii==1  
        p3=semilogy(4:0.1:8.1,10.^(Stat_log10N(:,4)),'b-','LineWidth',1.1);hold on
        p1=semilogy(mag_check,AnnualCumNum,'Color',[0.4 0.4 0.4],'LineWidth',1.5,'LineStyle','--','Marker','^','MarkerFaceColor','w'); hold on
        %plot(4:0.1:8.1,10.^(Stat_log10N(:,4)),'k-');   

    else
        p4=semilogy(4:0.1:8.1,10.^(Stat_log10N(:,4)),'k-','LineWidth',1.1);hold on
        p2=semilogy(mag_check,AnnualCumNum,'Color',[0.4 0.4 0.4],'LineWidth',1.5,'LineStyle','--','Marker','v','MarkerFaceColor','w'); hold on
        %plot(4:0.1:8.1,10.^(Stat_log10N(:,4)),'b-');
        %plot syncat_bg that is used in PSHA
        p5=semilogy(mag_check_syncat,bg_AnnualRate,'r-','LineWidth',1.1);hold on 
    end
       
end
 
set(gca,'fontsize',12)
legend([p1, p2, p3, p4, p5],{'Unfiltered SSA-GEM','Filtered SSA-GEM','Unfiltered fitted G-R','Filtered fitted G-R',['Rukwa-Malawi' newline 'Areal Source']}); grid on; 
    axis square; xlabel('Magnitude'); ylabel('Annual Frequency of Exceedance'); axis([Mc Mmax 0.001 10]); hold on;