%%%   Plot all GMPE responses for a hypothetical earthquake   %%%

clear
close all

mydir  = pwd; idcs   = strfind(mydir,'/');
addpath(mydir(1:idcs(end)-1));

load('syncat_PSHA_MSSD_input','GMPEindex','num_GMPEindex')

load GMPEcoef_Malawi

% GMPE options:
% 1)  Boore & Atkinson (2008) for shallow crustal events (NGA-WEST1)
% 2)  Campbell & Bozorgnia (2008) for shallow crustal events (NGA-WEST1)
% 3)  Chiou & Youngs (2008) for shallow crustal events (NGA-WEST1)
% 4)  Akkar & Bommer (2010) - Mw = 5.0-7.6; constant sigma model
% 5)  Boore, Stewart, Seyhan & Atkinson (2014) for shallow crustal events (NGA-WEST2)
% 6)  Akkar, Sandikkaya, & Bommer (2014) for shallow crustal events (Rjb model)
% 7)  Akkar, Sandikkaya, & Bommer (2014) for shallow crustal events (Repi model)
% 8)  Akkar, Sandikkaya, & Bommer (2014) for shallow crustal events (Rhypo model)
% 9)  Cauzzi, Faccioli, Vanini, & Bianchini (2015) for shallow crustal events
% 10) Chiou & Youngs (2014) for shallow crustal events (NGA-WEST2)
% 11) Atkinson and Adams (2013) for stable continental events (6th Seismic Hazard Model of Canada) - best/lower/upper branches 
% 12) Goulet et al. (2017) for stable continental events (NGA-East preliminary with site amplification; 6th Seismic Hazard Model of Canada) - models 1 to 13

%Specify event details to test GMM with

T = [0 0.1 0.2 0.3 0.5 1 2 3];
repi = 44.5;   vs30_site = 760; 

m=7.0; FMtype    = 0; Z10_Z25   = [0.1 0.8];  dip    = 53; fd=15;
num_T = length(T); rhypo     = sqrt(repi.^2 + fd.^2);  % Also used for Rrup
 
%% Plot GMPE 

GMmed=zeros(num_T,num_GMPEindex);

%NGA-East Fig Opt for if Goulet et al (2017) GMM selected (GMPE index 12)

%ngaeast_figopt=1 % median gm values from 10000 simulations through the logic tree
%ngaeast_gifopt=2 % weighted mean gm value for each weighted gmm branch

ngaeast_figopt=1;

figure(3);

for i=1:num_GMPEindex
    
     COEF = GMPEcoef_Malawi(GMPEindex(i),T);
    
    if GMPEindex(i)~=11 && GMPEindex(i)~=12
    
     %for active crust GMPE   
    [GMmed(:,i),Sintratmp,Sintertmp] = GMPE_Malawi(GMPEindex(i),COEF,m,repi,rhypo,0,repi,rhypo,vs30_site*ones(num_T,1),FMtype,fd,Z10_Z25,dip,T);     
     
    else %for nga-east2
        
        Tn=T; Vs=ones(1,num_T)*vs30_site;
        
       if ngaeast_figopt==1
        
       %assessed nga-east by running 10000 simulations through logic tree
       %and extracting mean value for each period
           
       GM_NumSimu=10000; 
        
       GMmed_tmp=zeros(GM_NumSimu,num_T);
        
        for ii = 1:length(Tn)
            
            for jj=1:GM_NumSimu
                
            lt_gmm = rand(1); 

            tn_pick = find(Tn(ii) == COEF.T);
        
            gmm_pick = find(lt_gmm<cumsum(COEF.LTweight(:,tn_pick)),1,'first');
    
            imtmp = COEF.table(:,2+tn_pick,:,gmm_pick);
            imtmp = reshape(imtmp(:),length(COEF.dist),length(COEF.mag),length(COEF.Vs30));
        
                if GMPEindex(i) == 11
                    immedtmp = interp3(COEF.maggrid,log10(COEF.distgrid),log10(COEF.vsgrid),imtmp,min(m,max(COEF.mag)),min(max(log10(rhypo),min(log10(COEF.dist))),max(log10(COEF.dist))),min(max(log10(vs30_site),min(log10(COEF.Vs30))),max(log10(COEF.Vs30))));
                elseif GMPEindex(i) == 12
                    immedtmp = interp3(COEF.maggrid,log10(COEF.distgrid),log10(COEF.vsgrid),imtmp,min(m,max(COEF.mag)),min(max(log10(repi),min(log10(COEF.dist))),max(log10(COEF.dist))),min(max(log10(vs30_site),min(log10(COEF.Vs30))),max(log10(COEF.Vs30))));
                end
                
            Sintra = sqrt((COEF.Sigma(tn_pick,1,gmm_pick).^2).*COEF.SigmaSplit(tn_pick,15)); % Natural logarithmic base
            Sinter = sqrt((COEF.Sigma(tn_pick,1,gmm_pick).^2).*COEF.SigmaSplit(tn_pick,14)); % Natural logarithmic base    
            
            GMmed_tmp(jj,ii)=((10.^immedtmp)/981)*exp(sqrt(Sintra.^2+Sinter.^2)*rand(1));
            
            
            end %end jj loop through GM_NumSimu
        
        GMmed(ii,i)=median(GMmed_tmp(:,ii));
        end %end ii loop through Tn
        
       else  %assessed nga-east through weighted mean of each GMM
           
           GMmed_tmp=zeros(size(COEF.LTweight,1),num_T);
           
           for ii=1:length(Tn)
               
                tn_pick = find(Tn(ii) == COEF.T);
                
                for ss=1:size(COEF.LTweight,1)
                    
                imtmp = COEF.table(:,2+tn_pick,:,ss);
                imtmp = reshape(imtmp(:),length(COEF.dist),length(COEF.mag),length(COEF.Vs30));
                
                if GMPEindex(i) == 11
                    immedtmp = interp3(COEF.maggrid,log10(COEF.distgrid),log10(COEF.vsgrid),imtmp,min(m,max(COEF.mag)),min(max(log10(rhypo),min(log10(COEF.dist))),max(log10(COEF.dist))),min(max(log10(vs30_site),min(log10(COEF.Vs30))),max(log10(COEF.Vs30))));
                elseif GMPEindex(i) == 12
                    immedtmp = interp3(COEF.maggrid,log10(COEF.distgrid),log10(COEF.vsgrid),imtmp,min(m,max(COEF.mag)),min(max(log10(repi),min(log10(COEF.dist))),max(log10(COEF.dist))),min(max(log10(vs30_site),min(log10(COEF.Vs30))),max(log10(COEF.Vs30))));
                end
                
                % take median gmm and multiply it by GMM branch logic tree
                % weighting
                GMmed_tmp(ss,ii)=(10.^immedtmp)/981*COEF.LTweight(ss,tn_pick);
                
               end %end ss loop  
               
               %sum weighted gmm to get weighted average
               GMmed(ii,i)=sum(GMmed_tmp(:,ii));
               
           end %end ii loop
      
       end
    end
    
    plot(T,GMmed(:,i),'Linewidth',1); hold on
    
    %clear COEF
end

xlabel('Period (s)'); ylabel('Ground Acceleration (g)')
legend({'Boore2014','Akkar2014','Chiou2014','Atkinson2013','Goulet2017'},'Location','northeast');
set(gca,'FontSize',13); axis square; grid on;
