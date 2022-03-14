%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%   Characteristic Earthquake Model of Youngs and Coppersmith (1985)   %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Youngs and Coppersmith (1985) 

% Formulations are based on Convertito et al. (2005).

% - Equation (3) -
% f(m) = 0 when m < Mmin and m > Mmax
% f(m) = beta*exp(-beta*(m-Mmin))/(1-exp(-beta*(Mmax-Mmin-DeltaM2)))*(1/(1+C)) for Mmin <= m <= Mc = Mmax - DeltaM2
% f(m) = beta*exp(-beta*(Mmax-Mmin-DeltaM1-DeltaM2))/(1-exp(-beta*(Mmax-Mmin-DeltaM2)))*(1/(1+C)) for Mc = Mmax - DeltaM2 <= m <= Mmax
% Mmin: Minimum magnitude
% Mmax: Maximum magnitude
% DeltaM1: 1.0 in Youngs and Coppersmith (1985)
% DeltaM2: 0.5 in Youngs and Coppersmith (1985)
% - Equation (4) -
% C = DeltaM2*beta*exp(-beta*(Mmax-Mmin-DeltaM1-DeltaM2))/(1-exp(-beta*(Mmax-Mmin-DeltaM2)))
% - Equation (5) -
% alphaC = alphaNC*DeltaM2*beta*exp(-beta*(Mmax-Mmin-DeltaM1-DeltaM2))/(1-exp(-beta*(Mmax-Mmin-DeltaM2))) = alphaNC*C
% - Equation (6) -
% alphaNC = rigidity*Af*S*(1-exp(-beta*(Mmax-Mmin-DeltaM2)))/(K*M0max*exp(-beta*(Mmax-Mmin-DeltaM2)))
% M0max = 10^(c*Mmax+d) = 10^(1.5*Mmax+9.05)
% - Equation (7) -
% K = b*10^(-c*DeltaM2)/(c-b) + b*exp(beta*DeltaM1)*(1-10^(-c*DeltaM2))/c
% - Equation (8) -
% alphaEXP = rigidity*Af*S*(c-b)*(1-exp(-beta*(Mmax-Mmin)))/(b*M0max*exp(-beta*(Mmax-Mmin)))

% Algorithm:
% - Specify: rigidity, Af, S, c, b, Mmax, Mmin, DeltaM1, and DeltaM2
% - Calculate alphaEXP using Equation (8)
% - Calculate alphaNC using Equations (6) and (7)
% - Calculate alphaC using Equation (5)
% - Calculate C using Equation (4)
% - Calculate f(m) using Equation (3)
% - Obtain F(m) and N(m)

% Input variables: para
% 1) Rigidity (Nm)
% 2) Fault plane area Af (m^2)
% 3) Slip rate S (m/year)
% 4) Gutenberg-Richter slope parameter b
% 5) Maximum magnitude Mmax
% 6) Minimum magnitude Mmin
% 7) Magnitude range for exponential distribution part DeltaM1
% 8) Magnitude range for characteristic distribution part DeltaM2

% Output variables:
% 1) alpha   - activity rates for non-characteristic events, characteristic events, and exponentially-distributed events
% 2) fm_char - probability density function, cumulative distribution function, and annual number of events exceeding m
% 3) fm_exp  - probability density function, cumulative distribution function, and annual number of events exceeding m

function [alpha,fm_char,fm_exp] = characteristic_magnitude_YC1985(para,mag_range,fignum,axis_range)

c     = 1.5;
beta  = log(10)*para(4);
M0max = 10^(c*para(5)+9.05);

K = para(4)*10^(-c*para(8))/(c-para(4)) + para(4)*exp(beta*para(7))*(1-10^(-c*para(8)))/c;

C = para(8)*beta*exp(-beta*(para(5)-para(6)-para(7)-para(8)))/(1-exp(-beta*(para(5)-para(6)-para(8))));

alpha(1) = para(1)*para(2)*para(3)*(1-exp(-beta*(para(5)-para(6)-para(8))))/(K*M0max*exp(-beta*(para(5)-para(6)-para(8))));

alpha(2) = alpha(1)*C;

alpha(3) = para(1)*para(2)*para(3)*(c-para(4))*(1-exp(-beta*(para(5)-para(6))))/(para(4)*M0max*exp(-beta*(para(5)-para(6))));

for ii = 1:length(mag_range)
    
    if mag_range(ii) < para(6) 
        fm_char(ii,1) = 0;
    elseif mag_range(ii) >= para(5)
        fm_char(ii,1) = 0;
    elseif para(6) <= mag_range(ii) && mag_range(ii) < para(5) - para(8)
        fm_char(ii,1) = beta*exp(-beta*(mag_range(ii)-para(6)))/(1-exp(-beta*(para(5)-para(6)-para(8))))*(1/(1+C));
    elseif para(5) - para(8) <= mag_range(ii) && mag_range(ii) < para(5)
        fm_char(ii,1) = beta*exp(-beta*(para(5)-para(6)-para(7)-para(8)))/(1-exp(-beta*(para(5)-para(6)-para(8))))*(1/(1+C));
    end
    
    if mag_range(ii) < para(6) 
        fm_exp(ii,1) = 0;
    elseif mag_range(ii) >= para(5)
        fm_exp(ii,1) = 0;
    else
        fm_exp(ii,1) = beta*exp(-beta*(mag_range(ii)-para(6)))/(1-exp(-beta*(para(5)-para(6))));
    end
    
end

fm_char(:,2) = min((mag_range(2)-mag_range(1))*cumsum(fm_char(:,1)),1);
fm_char(:,2) = fm_char(:,2)/fm_char(end,2);
fm_char(:,3) = (alpha(1)+alpha(2))*(1-fm_char(:,2));
fm_char(find(fm_char(:,1)==0),3) = NaN;

fm_exp(:,2) = min((mag_range(2)-mag_range(1))*cumsum(fm_exp(:,1)),1);
fm_exp(:,2) = fm_exp(:,2)/fm_exp(end,2);
fm_exp(:,3) = alpha(3)*(1-fm_exp(:,2));
fm_exp(find(fm_exp(:,1)==0),3) = NaN;

if ~isempty(fignum)
    
    disp(['Sum of area under PDF for the characteristic model = ',num2str(nanmax(fm_char(:,2)))]);
    disp(['Sum of area under PDF for the exponential model = ',num2str(nanmax(fm_exp(:,2)))]);
    
    figure (fignum); clf(fignum);
    subplot(221); plot(    mag_range,fm_char(:,1),'b-',mag_range,fm_exp(:,1),'r-'); hold on; axis square; title('PDF');
    subplot(222); plot(    mag_range,fm_char(:,2),'b-',mag_range,fm_exp(:,2),'r-'); hold on; axis square; title('CDF');
    subplot(212); semilogy(mag_range,fm_char(:,3),'b-',mag_range,fm_exp(:,3),'r-'); hold on; axis square; title('Recurrence rate'); axis(axis_range); 

end

