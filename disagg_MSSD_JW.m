%                                                        %
% M-function for implementing the deaggregation analysis %
%                                                        %

%MAX_PSA =AMAX_GM
%1) = GM
%2) = eps
%4) = mag
%8) = distance (rhypo)

function [] = disagg_MSSD(num_need,F,max_rad_dis,MAX_PSA,prob_level_deagg,Tn_pick,fignum,max_pmf,title_opt)

start = find(F >= 1-1/prob_level_deagg,1,'first'); %find index of events that occurred within prob level
%within AMAX_record

dmag = 0.2;
mag_bin = 4:dmag:8.4; %considers events from mag 4-8.5 at incremenets of 0.1
for ii = 1:length(mag_bin)-1
    mag_bin_center(ii) = (mag_bin(ii) + mag_bin(ii+1))/2;
end

ddis = 10;
dis_bin = 0:ddis:max_rad_dis; %max_rad set at 200 km, only considers events <200 km of site
for ii = 1:length(dis_bin)-1
    dis_bin_center(ii) = (dis_bin(ii) + dis_bin(ii+1))/2;
end

eps_bin        = [-inf -2.75 -2.25 -1.75 -1.25 -0.75 -0.25 0.25 0.75 1.25 1.75 2.25 2.75 inf];
eps_bin_center = [-3.0 -2.5 -2.0 -1.5 -1.0 -0.5 0.0 0.5 1.0 1.5 2.0 2.5 3.0];

transparancy = 0.7; 

COLORSCEME = colormap(parula(length(eps_bin)-1));
Bin_deagg = zeros(length(mag_bin)-1,length(dis_bin)-1,length(eps_bin)-1);

num_event_deagg = num_need - start + 1;%%Find how many events occurred for given probablilty 
% level within 1,000,000 long syn record, e.f. 2000 events for 1 in 500 year prob level)

disp(num2str(num_event_deagg))

for ii = 1:length(mag_bin)-1
    for jj = 1:length(dis_bin)-1
        for kk = 1:length(eps_bin)-1
            
            Bin_deagg(ii,jj,kk) = length(find( ...
                MAX_PSA(start:num_need,4) >= mag_bin(ii) & MAX_PSA(start:num_need,4) < mag_bin(ii+1) & ...
                MAX_PSA(start:num_need,3) >= dis_bin(jj) & MAX_PSA(start:num_need,3) < dis_bin(jj+1) & ...
                MAX_PSA(start:num_need,2) >= eps_bin(kk) & MAX_PSA(start:num_need,2) < eps_bin(kk+1)));
                %find number of events greater than each increment of mag 4-8, each increment of distance yp to 200 km, and
                %and epsilon (standard deviatiob variation within values
        end
    end           
end
        
%disp(['T = ',num2str(Tn_pick),' (s): ',num2str(num_event_deagg-sum(Bin_deagg(:))),' events are not included within the considered range']);
        
% Define the probability mass function.
PMF = Bin_deagg/num_event_deagg; %Number of events that fit criteria divded by total number of events
                    
% Mean magnitude, distance and epsilon
Mbar = 0;
for kk = 1:length(mag_bin)-1
    Mbar = Mbar + sum(sum(PMF(kk,:,:)))*mag_bin_center(kk);
end
        
Rbar = 0;
for kk = 1:length(dis_bin)-1
    Rbar = Rbar + sum(sum(PMF(:,kk,:)))*dis_bin_center(kk);
end
        
EPSbar = 0;
for kk = 1:length(eps_bin)-1
    EPSbar = EPSbar + sum(sum(PMF(:,:,kk)))*eps_bin_center(kk);
end
        
% Mode magnitude, distance and epsilon
Mode_index = zeros(4,1);
for ii = 1:length(mag_bin)-1
    for jj = 1:length(dis_bin)-1
        for kk = 1:length(eps_bin)-1
                
            if Mode_index(1) < PMF(ii,jj,kk)
                Mode_index(1) = PMF(ii,jj,kk);
                Mode_index(2) = ii;
                Mode_index(3) = jj;
                Mode_index(4) = kk;
            end
                    
        end
    end
end
        
Mmode   = mag_bin_center(Mode_index(2));
Rmode   = dis_bin_center(Mode_index(3));
EPSmode = eps_bin_center(Mode_index(4));

% Plot the M-D-epsilon deaggregation PMF.
figure(fignum)
Height = zeros(length(mag_bin)-1,length(dis_bin)-1);
for ii = 1:length(mag_bin)-1 % Magnitude
    for jj = 1:length(dis_bin)-1 % Distance
        for kk = 1:length(eps_bin)-1 % Epsilon

            if PMF(ii,jj,kk) > 0
                voxel([ dis_bin(1,jj)+ddis*0.1 mag_bin(1,ii)+dmag*0.1 Height(ii,jj)],[ddis*0.8 dmag*0.8 PMF(ii,jj,kk)],COLORSCEME(kk,:),transparancy); 
                hold on;
            end

            Height(ii,jj) = Height(ii,jj) + PMF(ii,jj,kk);
                
        end
    end
end
ylabel('Mw'); xlabel('Rupture distance (km)'); zlabel('PMF');
view(30,35); grid on; axis([min(dis_bin) max(dis_bin) min(mag_bin) max(mag_bin) 0 max_pmf]);
colormap(parula(length(eps_bin)-1));
set(gca,'fontsize',12);
set(get(gca,'xlabel'),'rotation',-13); set(get(gca,'ylabel'),'rotation',38); 
    
for ijk = 1:length(eps_bin)-1
    labels(ijk) = {num2str(eps_bin_center(ijk))};
end
lcolorbar(labels,'TitleString','\epsilon'); 
%lll.Label.String = '\epsilon';

%title revision by JW
title([title_opt(1) + " " + title_opt(2) + " " + title_opt(3) + " " + title_opt(5) + " g "  title_opt(4)],'fontweight','normal')
%title({['SA',num2str(Tn_pick),': ',num2str(round(1000*MAX_PSA(start,1))/1000),' g at ',num2str(prob_level_deagg),' RP years - # of data = ',num2str(num_event_deagg)]; ...
%    ['Mean(M,R,Eps)-',num2str(round(100*Mbar)/100),',',num2str(round(100*Rbar)/100),',',num2str(round(100*EPSbar)/100),' & Mode(M,R,Eps)-',num2str(Mmode),',',num2str(Rmode),',',num2str(EPSmode)]});
    
[Mbar Rbar EPSbar Mmode Rmode EPSmode]






