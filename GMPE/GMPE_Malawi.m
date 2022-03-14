%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%   Ground motion prediction equations for Malawi   %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%REVISED JW 23/06 SO GMMED CALCUATED AS MATRIX NOT IN FOR LOOP
%FOR GMPE OPTIONS 11&12. SEE GMPE_MALAWIARCHIVE FOR ORIGNAL FILE

% Output: 
% Median ground motion parameters are in g.
% Logarithmic standard deviation are natural logarithmic base.

% Sintra: within event variability
% Sinter: between event variability

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

function [GMmed,Sintra,Sinter] = GMPE_Malawi(GMPEindex,COEF,M,Rjb,Rrup,Rx,Repi,Rhypo,Vs,FMtype,Ztor,Z10_Z25,Dip,Tn)

Z10 = Z10_Z25(1)*1000; % m
Z25 = Z10_Z25(2); % km

num_GM = size(COEF,1);

if GMPEindex == 1
    
    % Mw: 5.0 - 8 and Rjb < 200 km
    % PGA and PSA (g) and PGV (cm/s) (cm) - geometric mean
    % ln(Y) = FM + FD + FS
    % FM = e1*U + e2*SS + e3*NS + e4*RS + e5*min(M-Mh,0) + e6*(min(M-Mh,0))^2 + e7*max(M-Mh,0)
    % FD = (c1 + c2*(M-Mref))*ln(sqrt(Rjb^2+h^2)/rref) + c3*(sqrt(Rjb^2+h^2)-rref)
    % FS = blin*ln(max(min(Vs30,1300),180)/760) + bnl*ln(0.06/0.1) + c*(ln(max(pga4nl,0.03)/0.03))^2 + d*(ln(max(pga4nl,0.03)/0.03))^3 for pga4nl <= 0.09
    % FS = blin*ln(max(min(Vs30,1300),180)/760) + bnl*ln(pga4nl/0.1) for pga4nl > 0.09
    % bnl = (b1-b2)*ln(min(max(Vs30,180),300)/300)/ln(180/300) + b2*ln(max(min(Vs30,760),300)/760)/ln(300/760) 
    % U: unspecified mechanism; SS: strike-slip mechanism; NS: normal-slip mechanism; RS: reverse-slip mechanism
    % Coef = [e1 e2 e3 e4 e5 e6 e7 Mh c1 c2 c3 h blin b1 b2 sigma tau-U sigma-TU tau-M sigma-TM]
    % Maximum T is 10.0 s        
    
    % "Unspecified" faulting type is not considered
    FM       = COEF(:,2); 
    FMpga4nl = -0.50350; 
    if FMtype == 1 % Normal
        FM       = COEF(:,3);
        FMpga4nl = -0.75472;
    elseif FMtype == 2 % Reverse
        FM       = COEF(:,4);
        FMpga4nl = -0.50970;
    end
    FM       = FM + COEF(:,5).*min(M-COEF(:,8),0) + COEF(:,6).*(min(M-COEF(:,8),0)).^2 + COEF(:,7).*max(M-COEF(:,8),0);
    FMpga4nl = FMpga4nl + 0.28805.*min(M-6.75,0) - 0.10164.*(min(M-6.75,0))^2;
    
    FD = (COEF(:,9)+COEF(:,10).*(M-4.5)).*log(sqrt(Rjb.^2+COEF(:,12).^2)/1) + COEF(:,11).*(sqrt(Rjb.^2+COEF(:,12).^2)-1);
    
    C = (3*log(0.09/0.06)-log(0.09/0.03))/(log(0.09/0.03))^2;
    D = -(2*log(0.09/0.06)-log(0.09/0.03))/(log(0.09/0.03))^3;
    for ii = 1:num_GM
        pga4nl(ii,1) = exp(FMpga4nl + (-0.66050 + 0.11970*(M-4.5))*log(sqrt(Rjb(ii,1).^2+1.350^2)/5) - 0.01151*(sqrt(Rjb(ii,1).^2+1.350^2)-5)); % Note rref = 5 is taken!
        bnl(ii,1)    = (COEF(ii,14)-COEF(ii,15))*log(min(max(Vs(ii,1),180),300)/300)/log(180/300) + COEF(ii,15)*log(max(min(Vs(ii,1),760),300)/760)/log(300/760);
        if pga4nl(ii,1) <= 0.09
            FS(ii,1) = COEF(ii,13)*log(max(min(Vs(ii,1),1300),180)/760) + bnl(ii,1).*(log(0.06/0.1) + C*(log(max(pga4nl(ii,1),0.03)/0.03)).^2 + D*(log(max(pga4nl(ii,1),0.03)/0.03)).^3);
        else
            FS(ii,1) = COEF(ii,13)*log(max(min(Vs(ii,1),1300),180)/760) + bnl(ii,1).*log(pga4nl(ii,1)/0.1);
        end
    end
    
    Sintra = COEF(:,16);
    Sinter = COEF(:,19);
        
    GMmed = exp(FM + FD + FS); % (g)   
        
elseif GMPEindex == 2
    
    % Mw: from 4.0 to 7.5-8.5 and Rrup < 200 km
    % PGA and PSA (g), PGV (cm/s), and PGD (cm) - geometric mean
    % ln(Y) = fmag + fdis + fflt + fhng + fsite + fsed
    % fmag = c0 + c1*M + c2*max((M-5.5),0) + c3*max((M-6.5),0)
    % fdis = (c4+c5*M)*ln(sqrt(Rrup^2+c6^2))
    % fflt = c7*Frv*min(Ztor,1) + c8*Fnm
    % fhng = c9*fhngr*min(max(2*(M-6),0),1)*max((20-Ztor)/20,0)*min((90-delta)/20,1)
    % fsite = c10*ln(Vs30/k1) + k2*(ln(A1100+1.88*(Vs30/k1)^1.18)-ln(A1100+1.88)) for Vs30 < k1
    % fsite = (c10+k2*1.18)*ln(min(Vs30,1100)/k1) for Vs30 >= k1
    % fsed = c11*(Z25-1) for Z25 < 1; 0 for 1 <= Z25 <= 3; c12*k3*exp(-0.75)*(1-exp(-0.25*(Z25-3))) for Z25 > 3
    % Frv = 1 for reverse and reverse-oblique events (30 < rake angle < 150)
    % Fnm = 1 for normal and normal-oblique events (-150 < rake angle < -30)
    % Ztor = depth to the top of the coseismic rupture plane (km)
    % A1100 = median estimate of PGA on a reference rock outcrop.
    % Z25 = depth to the 2.5 km/s shear-wave velocity horizon
    % sigma_total = sqrt(sigma^2 + tau_lnY^2) (note: sigma_c may be included)
    % sigma = sqrt(sigma_lnY^2 + alpha^2*(sigma_lnPGA^2 - sigma_lnAF^2) + 2*alpha*rho*sqrt(sigma_lnY^2-sigma_lnAF^2)*sqrt(sigma_lnPGA^2 - sigma_lnAF^2))
    % alpha = k2*A1100*(1/(A1100+1.88*(min(Vs30,k1)/k1)^1.18) - 1/(A1100+1.88))
    % Coef = [c0 c1 c2 c3 c4 c5 c6 c7 c8 c9 c10 c11 c12 k1 k2 k3 sigma-lnY tau-lnY sigma-c sigma-T sigma-Arb rho]
    % Maximum T is 10.0 s    
        
    FMAG = COEF(:,1) + COEF(:,2)*M + COEF(:,3)*max((M-5.5),0) + COEF(:,4)*max((M-6.5),0);
    
    FDIS = (COEF(:,5)+COEF(:,6)*M).*log(sqrt(Rrup.^2+COEF(:,7).^2));
    
    FFLT      = 0;
    FFLTA1100 = 0;
    if FMtype == 1 % Normal
        FFLT      = COEF(:,9);
        FFLTA1100 = -0.120;
    elseif FMtype == 2 % Reverse
        FFLT      = COEF(:,8)*min(Ztor,1);
        FFLTA1100 = 0.280*min(Ztor,1);
    end
    
    tmp = find(Rjb > 0);
    fhngr = ones(num_GM,1);
    if Ztor < 1
        fhngr(tmp,1) = (max(Rrup(tmp,1),sqrt(Rjb(tmp,1).^2+1))-Rjb(tmp,1))./max(Rrup(tmp,1),sqrt(Rjb(tmp,1).^2+1));
    elseif Ztor >= 1
        fhngr(tmp,1) = (Rrup(tmp,1)-Rjb(tmp,1))./Rrup(tmp,1);
    end
    FHNG = COEF(:,10).*fhngr*min(max(2*(M-6),0),1)*max((20-Ztor)/20,0)*min((90-Dip)/20,1);
    
    FSED      = 0;
    FSEDA1100 = 0;
    if Z25 < 1
        FSED      = COEF(:,12)*(Z25-1);
        FSEDA1100 = 0.040*(Z25-1);
    elseif Z25 > 3
        FSED      = COEF(:,13)*COEF(:,16)*exp(-0.75)*(1-exp(-0.25*(Z25-3)));
        FSEDA1100 = 0.610*1.839*exp(-0.75)*(1-exp(-0.25*(Z25-3)));
    end

    A1100 = exp(-1.715+0.5*M-0.53*max((M-5.5),0)-0.262*max((M-6.5),0) + (-2.118+0.17*M)*log(sqrt(Rrup.^2+5.6^2)) + FFLTA1100 + ...
        0.490*fhngr*min(max(2*(M-6),0),1)*max((20-Ztor)/20,0)*min((90-Dip)/20,1) + FSEDA1100 + (1.058-1.186*1.18)*log(1100/865));
    if length(A1100) < num_GM; A1100 = A1100*ones(num_GM,1); end
    sigma_lnAF = 0.3;
    for ii = 1:num_GM
        
        if Vs(ii,1) < COEF(ii,14)
            FSITE(ii,1) = COEF(ii,11)*log(Vs(ii,1)/COEF(ii,14)) + COEF(ii,15)*(log(A1100(ii)+1.88*(Vs(ii,1)/COEF(ii,14))^1.18)-log(A1100(ii)+1.88));
        else
            FSITE(ii,1) = (COEF(ii,11)+COEF(ii,15)*1.18)*log(min(Vs(ii,1),1100)/COEF(ii,14));
        end
    
    end
    
    % Orientation variability is not included.
    Alpha = A1100.*COEF(:,15).*(1./(A1100+1.88*(min(Vs,COEF(:,14))./COEF(:,14)).^1.18) - 1./(A1100+1.88));
    
    Sintra = sqrt(COEF(:,17).^2 + Alpha.^2*(0.478^2 - sigma_lnAF^2) + 2*Alpha.*COEF(:,22).*sqrt(COEF(:,17).^2-sigma_lnAF^2)*sqrt(0.478^2 - sigma_lnAF^2));
    Sinter = COEF(:,18);

    GMmed = exp(FMAG + FDIS + FFLT + FHNG + FSED + FSITE); % (g)        
    
elseif GMPEindex == 3
    
    % Mw: 4.0 - 8.0-8.5 and Rrup < 200 km
    % PGA and PSA (g) and PGV (cm/s) - geometric mean
    % ln(y) = ln(yref) + phi1*min(ln(Vs30/1130),0) + phi2*(exp(phi3*(min(Vs30,1130)-360)) - exp(phi3*770))*ln((yref+phi4)/phi4) + phi5*(1-1/cosh(phi6*max(Z10-phi7,0))) + phi8/cosh(0.15*max(Z10-15,0))
    % ln(yref) = c1 + (c1a*Frv + c1b*Fnm + c7*(Ztor-4))*(1-AS) + (c10 + c7a*(Ztor-4))*AS + c2*(M-6) + ((c2-c3)/cn)*ln(1+exp(cn*(cm-M))) + c4*ln(Rrup + c5*cosh(c6*max(M-chm,0))) + (c4a-c4)*ln(sqrt(Rrup^2+crb^2)) + (cg1 + cg2/cosh(max(M-cg3,0)))*Rrup + c9*Fhw*tanh(Rx*cos(delta)*cos(delta)/c9a)*(1-sqrt(Rjb^2+Ztor^2)/(Rrup+0.001))
    % sigma-T^2 = (1+NL0)^2*tau^2 + sigma-nl0^2
    % NL0 = phi2*(exp(phi3*(min(Vs30,1130)-360)) - exp(phi3*770))*yref/(yref+phi4)
    % sigma-nl0 = (sigma1 + 0.5*(sigma2-sigma1)*(min(max(M,5),7)-5) + sigma4*AS)*sqrt(sigma3*Finfer + 0.7*Fmeasure + (1+NL0)^2)
    % tau = tau1 + 0.5*(tau2-tau1)*(min(max(M,5),7)-5)
    % Rx (km): site coordinate measured perpendicular to the fault strike from the surface projection of the updip edge of the fault rupture
    % Fhw: hanging-wall flag - "1" for Rx >= 0 and "0" for Rx < 0 
    % AS: aftershock flag
    % Z10 (m): depth to shear wave velocity of 1.0 km/s
    % Finfer = 1 if Vs30 is inferred from geology, and Fmeasure = 1 if Vs30 is measured
    % Coef = [c1 c1a c1b cn cm c5 c6 c7 c7a c9 c9a c10 cg1 cg2 phi1 phi2 phi3 phi4 phi5 phi6 phi7 phi8 tau1 tau2 sigma1 sigma2 sigma3 sigma4]
    % Additional constants: c2 = 1.06, c3 = 3.45, c4 = -2.1, c4a = -0.5, crb = 50, chm = 3, cg3 = 4
    % Maximum T is 10.0 s   
        
    AS       = 0; 
    Finfer   = 1; 
    Fmeasure = 0;
    c3       = 3.45; 
    
    for ii = 1:num_GM
        if Rx(ii,1) >= 0
            Fhw(ii,1) = 1;
        else
            Fhw(ii,1) = 0;
        end
    end

    Fnm = 0; 
    Frv = 0;
    if FMtype == 1 % Normal
        Fnm = 1;
    elseif FMtype == 2 % Reverse
        Frv = 1;
    end

    % Coef = [c1 c1a c1b cn cm c5 c6 c7 c7a c9 c9a c10 cg1 cg2 phi1 phi2 phi3 phi4 phi5 phi6 phi7 phi8 tau1 tau2 sigma1 sigma2 sigma3 sigma4]
    % Additional constants: c2 = 1.06, c3 = 3.45, c4 = -2.1, c4a = -0.5, crb = 50, chm = 3, cg3 = 4
    lnYref = COEF(:,1) + (Frv*COEF(:,2) + Fnm*COEF(:,3) + COEF(:,8)*(Ztor-4))*(1-AS) + (COEF(:,12) + COEF(:,9)*(Ztor-4))*AS + ...
        1.06*(M-6) + ((1.06-c3)./COEF(:,4)).*log(1+exp(COEF(:,4).*(COEF(:,5)-M))) - 2.1*log(Rrup + COEF(:,6).*cosh(COEF(:,7)*max(M-3,0))) + ...
        1.6*log(sqrt(Rrup.^2+50^2)) + (COEF(:,13) + COEF(:,14)/cosh(max(M-4,0))).*Rrup + ...
        Fhw.*COEF(:,10).*tanh(Rx*cos(Dip/180*pi)*cos(Dip/180*pi)./COEF(:,11)).*max((1-sqrt(Rjb.^2+Ztor^2)./(Rrup+0.001)),0);
    
    GMmed = exp(lnYref + COEF(:,15).*min(log(Vs/1130),0) + COEF(:,16).*(exp(COEF(:,17).*(min(Vs,1130)-360)) - exp(COEF(:,17)*770)).*log((exp(lnYref)+COEF(:,18))./COEF(:,18)) + ...
        COEF(:,19).*(1-1./cosh(COEF(:,20).*max(Z10-COEF(:,21),0))) + COEF(:,22)/cosh(0.15*max(Z10-15,0))); % (g)
    
    NL0 = COEF(:,16).*(exp(COEF(:,17).*(min(Vs,1130)-360)) - exp(COEF(:,17)*770)).*exp(lnYref)./(exp(lnYref)+COEF(:,18));
    
    Sintra = (COEF(:,25) + 0.5*(COEF(:,26)-COEF(:,25))*(min(max(M,5),7)-5) + COEF(:,28)*AS).*sqrt(COEF(:,27)*Finfer + 0.7*Fmeasure + (1+NL0).^2);%eq 20 within event variance
    Sinter = (1+NL0).*(COEF(:,23) + 0.5*(COEF(:,24)-COEF(:,23))*(min(max(M,5),7)-5));% eq 19, between event variance
    
elseif GMPEindex == 4
    
    % Akkar and Bommer (2010) - Mw = 5.0-7.6
    % log10(PGV/PGA/SD) = b1 + b2*M + b3*M^2 + (b4+b5*M)*log10(sqrt(Rjb^2+b6^2)) + b7*Ss + b8*Sa + b9*Fn + b10*Fr
    % Ss = 1: soft soil site index (Vs30 < 360 m/s)
    % Sa = 1: firm soil site index (360 < Vs30 < 750 m/s)
    % Fn = 1: normal faulting index
    % Fr = 1: reverse faulting index
    % Constant sigma model
    % Coef = [b1 b2 b3 b4 b5 b6 b7 b8 b9 b10 sigma-intra sigma-inter sigma-total]
    
    for ii = 1:num_GM
        Soil(ii,1) = 0;
        if Vs(ii,1) < 360
            Soil(ii,1) = COEF(ii,7);
        elseif Vs(ii,1) >= 360 && Vs(ii,1) <= 750
            Soil(ii,1) = COEF(ii,8);
        end
    end
    
    Fault = 0; 
    if FMtype == 1 % Normal
        Fault = COEF(:,9);
    elseif FMtype == 2 % Reverse
        Fault = COEF(:,10);
    end
    
    GMmed = 10.^(COEF(:,1) + COEF(:,2)*M + COEF(:,3)*M*M + (COEF(:,4)+COEF(:,5)*M).*log10(sqrt(Rjb.^2+COEF(:,6).^2)) + Soil + Fault)/981; % (g)    

    Sintra = COEF(:,11)*log(10); % Natural logarithmic base
    Sinter = COEF(:,12)*log(10); % Natural logarithmic base
    
elseif GMPEindex == 5

    % Orientation-independent measure (no conversion) 
    % log(Y) = FE + FP + FS + epsilon
    % FE = e0*U + e1*SS + e2*NS + e3*RS + e4*(M-Mh) + e5*(M-Mh)^2 (M <= Mh)
    % FE = e0*U + e1*SS + e2*NS + e3*RS + e6*(M-Mh) (M > Mh)
    % FP = (c1 + c2*(M-4.5))*log(R) + c3*(R-1) (global model only)
    % R = sqrt(Rjb^2 + h^2)
    % FS = lnFlin + lnFnl (no basin effect)
    % lnFlin = c*log(min(Vs30,Vc)/760)
    % lnFnl  = [f4*(exp(f5*[min(Vs30,760)-360]) - exp(f5*400))]*log((PGAr+0.1)/0.1)
    % sigma = sqrt(sigma_intra^2 + sigma_inter^2)
    % sigma_inter = tau1 + (tau2-tau1)*min(max(0,M-4.5),1)
    % sigma_intra = phi1 + (phi2-phi1)*min(max(0,M-4.5),1) + DphiR*min(max([log(Rjb/R1)/log(R2/R1)],0),1) - DphiV*min(max([log(V2/Vs30)/log(V2/V1)],0),1) 
    % Y is PGA/PSA (g) or PGV (cm/s)
    % U: unspecified mechanism; SS: strike-slip mechanism; NS: normal-slip mechanism; RS: reverse-slip mechanism
    % Coef = [e0 e1 e2 e3 e4 e5 e6 Mh c1 c2 c3 h Dc3G Dc3CT Dc3JI c Vc f4 f5 f6 f7 R1 R2 DfR DfV V1 V2 f1 f2 t1 t2]
    % Maximum T is 10.0 s        
    
    % Sigma should be evaluated using Rjb
    Sintra = COEF(:,28) + (COEF(:,29)-COEF(:,28))*min(max(0,M-4.5),1) + ...
        COEF(:,24).*min(max((log(Rjb./COEF(:,22))./log(COEF(:,23)./COEF(:,22))),0),1) - ...
        COEF(:,25).*min(max((log(COEF(:,27)./Vs)./log(COEF(:,27)./COEF(:,26))),0),1);
    Sinter = COEF(:,30) + (COEF(:,31)-COEF(:,30))*min(max(0,M-4.5),1);
    
    % PGAr
    FEpgar = 0.4856; %for strike-slip mechanism
    if FMtype == 1 % Normal
        FEpgar = 0.2459;
    elseif FMtype == 2 % Reverse
        FEpgar = 0.4539;
    end
    if M <= 5.5
        FEpgar = FEpgar + 1.431*(M-5.5) + 0.05053*(M-5.5)^2; 
    else
        FEpgar = FEpgar - 0.1662*(M-5.5);
    end
    FPpgar = (-1.134 + 0.1917*(M-4.5))*log(sqrt(Rjb.^2+4.5^2)) - 0.008088*(Rjb-1);
    PGAr = exp(FEpgar + FPpgar);

    % "Unspecified" faulting type is not considered!
    FE = COEF(:,2); 
    if FMtype == 1 % Normal
        FE = COEF(:,3);
    elseif FMtype == 2 % Reverse
        FE = COEF(:,4);
    end
   
    % Magnitude-dependent term
    for ii = 1:num_GM
        if M <= COEF(ii,8)
            FE(ii,1) = FE(ii,1) + COEF(ii,5)*(M-COEF(ii,8)) + COEF(ii,6)*(M-COEF(ii,8))^2;
        else
            FE(ii,1) = FE(ii,1) + COEF(ii,7)*(M-COEF(ii,8));
        end
    end
        
    R = sqrt(Rjb.^2 + COEF(:,12).^2);
    
    % Distance-dependent term (no regional dependency)
    FP = (COEF(:,9)+COEF(:,10)*(M-4.5)).*log(R) + COEF(:,11).*(R-1);
    
    % Site-dependent term
    for ii = 1:num_GM
        Vs2(ii,1) = min(Vs(ii,1),COEF(ii,17)); 
    end
    FS = COEF(:,16).*log(Vs2/760) + (COEF(:,18).*(exp(COEF(:,19).*(min(Vs,760)-360)) - exp(COEF(:,19)*400))).*log((PGAr+0.1)/0.1);

    GMmed = exp(FE + FP + FS); % (g)
            
elseif GMPEindex == 6 || GMPEindex == 7 || GMPEindex == 8
    
    % ln(Y) = ln(Yref(Mw,R,Fn,Fr)) + ln(S(Vs,PGAref)) + epsilon*sigma	
    % ln(Yref) = a1 + 0.0029*(Mw-6.75) + a3*(8.5-Mw)^2 + (a4+0.2529*(Mw-6.75))*ln(sqrt(R^2+7.5^2)) + a8*Fn + a9*Fr   for Mw <= 6.75	
    % ln(Yref) = a1 - 0.5096*(Mw-6.75) + a3*(8.5-Mw)^2 + (a4+0.2529*(Mw-6.75))*ln(sqrt(R^2+7.5^2)) + a8*Fn + a9*Fr   for Mw > 6.75	
    % ln(S) = b1*ln(Vs/750) + b2*ln([PGAref+2.5*(Vs/750)^3.2]/[(PGAref+2.5)*(Vs/750)^3.2])   for Vs <= 750	
    % ln(S) = b1*ln(min(Vs,1000)/750) for Vs > 750	
    % PGA and PSA (g) and PGV (cm/s)	
    % R can be Rjb, Repi, Rhypo	
    % Y is geometric mean	
    % Mw: 4.0 to 7.6 and Rjb: 0 to 200 km	
    % Coef = [a1 a3 a4 a8 a9 b1 b2 sigma1 sigma2 sigma]
    % T: 0 to 4.0 s	
    % For ss mechanisms, Fnm and Frv should equal 0
    
    Fnm = 0; 
    Frv = 0;
    if FMtype == 1 % Normal
        Fnm = 1;
    elseif FMtype == 2 % Reverse
        Frv = 1;
    end
    
    if GMPEindex == 6
        R = Rjb;
        if M <= 6.75
            PGAr = exp(1.85329 + 0.0029*(M-6.75) - 0.02807*(8.5-M)^2 + (-1.23452+0.2529*(M-6.75))*log(sqrt(R.^2+7.5^2)) - 0.1091*Fnm + 0.0937*Frv);
        else
            PGAr = exp(1.85329 - 0.5096*(M-6.75) - 0.02807*(8.5-M)^2 + (-1.23452+0.2529*(M-6.75))*log(sqrt(R.^2+7.5^2)) - 0.1091*Fnm + 0.0937*Frv);
        end
    elseif GMPEindex == 7
        R = Repi;
        if M <= 6.75
            PGAr = exp(2.52977 + 0.0029*(M-6.75) - 0.05496*(8.5-M)^2 + (-1.31001+0.2529*(M-6.75))*log(sqrt(R.^2+7.5^2)) - 0.1091*Fnm + 0.0937*Frv);
        else
            PGAr = exp(2.52977 - 0.5096*(M-6.75) - 0.05496*(8.5-M)^2 + (-1.31001+0.2529*(M-6.75))*log(sqrt(R.^2+7.5^2)) - 0.1091*Fnm + 0.0937*Frv);
        end
    elseif GMPEindex == 8
        R = Rhypo;
        if M <= 6.75
            PGAr = exp(3.26685 + 0.0029*(M-6.75) - 0.04846*(8.5-M)^2 + (0.2529+0.2529*(M-6.75))*log(sqrt(R.^2+7.5^2)) - 0.1091*Fnm + 0.0937*Frv);
        else
            PGAr = exp(3.26685 - 0.5096*(M-6.75) - 0.04846*(8.5-M)^2 + (0.2529+0.2529*(M-6.75))*log(sqrt(R.^2+7.5^2)) - 0.1091*Fnm + 0.0937*Frv);
        end
    end

    if M <= 6.75
        FEP = COEF(:,1) + 0.0029*(M-6.75) + COEF(:,2)*(8.5-M)^2 + (COEF(:,3)+0.2529*(M-6.75)).*log(sqrt(R.^2+7.5^2)) + COEF(:,4)*Fnm + COEF(:,5)*Frv;
    else
        FEP = COEF(:,1) - 0.5096*(M-6.75) + COEF(:,2)*(8.5-M)^2 + (COEF(:,3)+0.2529*(M-6.75)).*log(sqrt(R.^2+7.5^2)) + COEF(:,4)*Fnm + COEF(:,5)*Frv;
    end
        
    % Site-dependent term
    b2 = zeros(num_GM,1);
    sitesoft = find(Vs<=750);
    b2(sitesoft) = COEF(sitesoft,7);
    FS = COEF(:,6).*log(min(Vs,1000)/750) + b2.*log((PGAr+2.5*(Vs/750).^3.2)./((PGAr+2.5).*(Vs/750).^3.2));

    GMmed = exp(FEP + FS); % (g)
    
    Sintra = COEF(:,8);
    Sinter = COEF(:,9);
    
elseif GMPEindex == 9
    
    % log10(Y) = fM + fR + fS + fFM + epsilon	
    % fM = c1 + m1*Mw + m2*Mw^2	
    % fR = (r1 + r2*Mw)*log10(Rrup + r3)	
    % fS = sB*SB + sC*SC + sD*SD	
    % fS = bv*log10(Vs/VA)	
    % fS = bv800*log10(Vs/800)	
    % fFM = fN*FN + fR*FR + fSS*FSS	
    % Y is 5%-damped displacement response spectrum in cm. It requires a conversion to obtain psuedo-spectral acceleration. PGA (cm/s/s)	
    % Y is geometric mean	
    % Mw: 4.5 to 7.9 and Rrup: 0 to 150 km	
    % Coef = [c1 m1 m2 r1 r2 r3 sB sC sD bV bV800 VA fN fR fSS phi tau sigma tauM sigmaM]
    % T: 0 to 10.0 s	
        
    Fss = 0; 
    Fnm = 0; 
    Frv = 0;
    if FMtype == 0 % Strike-slip
        Fss = 1;
    elseif FMtype == 1 % Normal
        Fnm = 1;
    elseif FMtype == 2 % Reverse
        Frv = 1;
    end
    
    FM = COEF(:,1) + COEF(:,2)*M + COEF(:,3)*M*M + Fnm*COEF(:,13) + Frv*COEF(:,14) + Fss*COEF(:,15);
    
    FR = (COEF(:,4)+COEF(:,5)*M).*log10(Rrup+COEF(:,6));
    
    FSoption = 3; % 1) Discrete site class, 2) Variable reference Vs, 3) Constant reference Vs of 800 m/s
    if FSoption == 1
        FS = zeros(num_GM,1);
        vsclassB = find(Vs >= 360 & Vs < 800);
        vsclassC = find(Vs >= 180 & Vs < 360);
        vsclassD = find(Vs < 180);
        FS(vsclassB) = COEF(vsclassB,7);
        FS(vsclassC) = COEF(vsclassC,8);
        FS(vsclassD) = COEF(vsclassD,9);
    elseif FSoption == 2
        FS = COEF(:,10).*log10(Vs./COEF(:,12));
    elseif FSoption == 3
        FS = COEF(:,11).*log10(Vs/800);
    end

    GMmed = 10.^(FM + FR + FS); % Spectral displacement (cm)
    
    Tn(find(Tn == 0)) = 0.01;
    
    GMmed = (GMmed.*(2*pi./Tn).^2)/981; % Pseudo-spectral acceleration (g)
    
    Sintra = COEF(:,16)*log(10); % Natural logarithmic base
    Sinter = COEF(:,19)*log(10); % Natural logarithmic base
    
elseif GMPEindex == 10
    
    % Rx (km): site coordinate measured perpendicular to the fault strike from the surface projection of the updip edge of the fault rupture
    % Fhw: hanging-wall flag - "1" for Rx >= 0 and "0" for Rx < 0 
    % GMmed defined by equation 12 in C&Y14
    % Sinter & Sintra defined by equation 13 in C&Y14
    % Sinter = Tau_1+(Tau_2-Tau1)/1.5(min(max(M,5),6.5)-5)
    % Sintra = sigma1 + (sigma2-sigma1)*min(max(M,5),6.5)-5) * (sigma3*F_inferred+0.7*F_measured+(1+NL_0)^2)^(0.5)

    if Rx >= 0
        Fhw = 1;
    else
        Fhw = 0;
    end
   
    Fnm = 0; Frv = 0;
    if FMtype == 1 % Normal
        Fnm = 1;
    elseif FMtype == 2 % Reverse
        Frv = 1;
    end
    
    % Depth to 1 km/s as a function of vs30 (eq. 1 in C&Y14)
    z_1 = exp(-7.15/4*log((Vs.^4+570.94^4)/(1360^4+570.94^4)));
    
    if Z10 ==999
        d_Z1 = 0;
    else
        d_Z1 = Z10 - z_1; % note z10 previously converted to m
    end
      
    term1 = COEF(:,2); % Parameter c1
        
    % Style of faulting term for c1_a to c1_d
    term2(:,1) = (COEF(:,3)+COEF(:,5)/(cosh(2*max(M-4.5,0))))*Frv;
    term3(:,1) = (COEF(:,4)+COEF(:,6)/cosh(2*max(M-4.5,0)))*Fnm;
        
    % Ztor term (i.e. depth to top of rupture plane)
    if Frv == 1
        E_Ztor = (max(2.704-1.226*max(M-5.849,0),0))^2;% eq 4 in C&Y14
    else
        E_Ztor = (max(2.673-1.136*max(M-4.970,0),0))^2;% eq 5 in C&Y14
    end
    
    if Ztor == 999
        Ztor = E_Ztor;
    end
    delta_ZTOR = Ztor - E_Ztor;
    
    % use coef c7 and c7_b in C&Y14
    term4(:,:) = (COEF(:,7)+COEF(:,8)/cosh(2*max(M-4.5,0))).*delta_ZTOR;
      
    % Dip term using coef c11 and c11_b in C&Y14
    term5(:,:) = (COEF(:,9)+COEF(:,10)/cosh(2*max(M-4.5,0)))*(cos(Dip*pi()/180.0)^2);
    
    % fmag using coef c2, c3, cn, and cm
    term6(:,:) = COEF(:,11)*(M-6);
    term7(:,:) = (COEF(:,11)-COEF(:,12))./COEF(:,13).*log(1+exp(COEF(:,13).*(COEF(:,14)-M)));
        
    % Distance Scaling and attenuation term
    % Incorporate coef c4, c5, c6, and c_HM
    term8(:,:) = COEF(:,15).*log(Rrup+COEF(:,16).*cosh(COEF(:,17).*max(M-COEF(:,18),0)));
    % Incorporate coef c4, c4a, and c_RB
    term9(:,:) = (COEF(:,19)-COEF(:,15)).*log(sqrt(Rrup.^2+COEF(:,20).^2));
    % Inorporate coef cg_1 to cg_3
    term10(:,:) = (COEF(:,21)+COEF(:,22)./(cosh(max(M-COEF(:,23),0)))).*Rrup;
      
    % Directivity
    d_DPP = 0; % for median estimate
    term11(:,:) = COEF(:,24).*max(1-max(Rrup-40,0)/30,0)*min(max(M-5.5,0)/0.8,1).*exp(-COEF(:,25).*(M-COEF(:,26)).^2).*d_DPP;
    
    % Hanging wall term
    term12(:,:) = COEF(:,27).*Fhw.*cos(Dip*pi()/180.0).*(COEF(:,28)+(1-COEF(:,28)).*tanh(Rx./COEF(:,29))).*(1-sqrt(Rjb.^2+Ztor.^2)./(Rrup+1));
        
    ln_yrefij(:,:) = term1 + term2 + term3 + term4 + term5 + term6 + term7 + term8 + term9 + term10 + term11 + term12;
        
    yrefij(:,:)=exp(ln_yrefij);
        
    % Site response % coef phi1 to phi6
    term14(:,:) = COEF(:,30).*min(log(Vs/1130),0);
    term15(:,:) = COEF(:,31).*(exp(COEF(:,32).*(min(Vs,1130)-360))-exp(COEF(:,32).*(1130-360))).*log((yrefij+COEF(:,33))./COEF(:,33));
    term16(:,:) = COEF(:,34).*(1-exp(-d_Z1./COEF(:,35)));

    GMmed = yrefij.*exp(term14 + term15 + term16);
        
    % Compute standard deviation % vs30 option is inferred.
    % Use eq 13 in CY14  
    NL0 = COEF(:,31).*(exp(COEF(:,32).*(min(Vs,1130)-360))-exp(COEF(:,32).*(1130-360))).*(yrefij./(yrefij+COEF(:,33)));
    Sintra = (COEF(:,36) + (COEF(:,37) - COEF(:,36))./1.5.*(min(max(M,5),6.5)-5)).*sqrt((COEF(:,38).*1) + (1+NL0).^2); % sigmaNL_0 in CY14 (within event variance) 
    sinter_tmp = COEF(:,39) + (COEF(:,40) - COEF(:,39))/1.5.*(min(max(M,5),6.5)-5);% Tau in CY14 (between event variability)
    % sigma = sqrt((1+NL0).^2.*Sintra.^2 + Sintra.^2); % total event variance

    Sinter = (1+NL0).*sinter_tmp; %addition from Katsu, see email 16/06/21 
    
elseif GMPEindex == 11 || GMPEindex == 12
    
    % COEF = CanadaSHM6_2015/CanadaSHM6_2020
    % This structure variable includes: table, T, mag, dist, Vs30, Sigma, SigmaSplit, LTweight, maggrid, distgrid, vsgrid
    % The distance measure is hypocentral distance for CanadaSHM6_2015.
    % The distance measure is rupture distance for CanadaSHM6_2020.
    % SigmaSplit contains the 'average' variance ratios of inter-event variability (column 14) and intra-event variability (column 15).
    % The average of the variance ratios is computed based on three values for M5, M6, and M7+ (see Table 5.5. of Goulet et al. (2017)) - this is consistent with the sigma model adopted in CanadaSHM6 
    
    lt_gmm = rand(1);
    %lt_gmm = 0.5;% FIX THIS VALUE FOR TESTING
    
    tn_pick=zeros(length(Tn),1); gmm_pick=zeros(length(Tn),1); immedtmp = zeros(length(Rhypo),length(Tn)); 
    Sintra = zeros(length(Rhypo),length(Tn)); Sinter = zeros(length(Rhypo),length(Tn)); 
    
    M=ones(length(Rhypo),1)*M; %repeat M. so every array in Interp3 function has size (4,1)
    
    for ii = 1:length(Tn)
        
        tn_pick(ii,:) = find(Tn(ii) == COEF.T);
        
        gmm_pick(ii,:) = find(lt_gmm<cumsum(COEF.LTweight(:,tn_pick(ii,:))),1,'first');
    
        imtmp = COEF.table(:,2+tn_pick(ii,:),:,gmm_pick(ii,:));
        imtmp = reshape(imtmp(:),length(COEF.dist),length(COEF.mag),length(COEF.Vs30));
        if GMPEindex == 11
            immedtmp(:,ii) = interp3(COEF.maggrid,log10(COEF.distgrid),log10(COEF.vsgrid),imtmp,min(M,max(COEF.mag)),min(max(log10(Rhypo),min(log10(COEF.dist))),max(log10(COEF.dist))),min(max(log10(Vs),min(log10(COEF.Vs30))),max(log10(COEF.Vs30))));
        elseif GMPEindex == 12
            immedtmp(:,ii) = interp3(COEF.maggrid,log10(COEF.distgrid),log10(COEF.vsgrid),imtmp,min(M,max(COEF.mag)),min(max(log10(Rrup),min(log10(COEF.dist))),max(log10(COEF.dist))),min(max(log10(Vs),min(log10(COEF.Vs30))),max(log10(COEF.Vs30))));
        end
        
        Sintra(:,ii) = sqrt((COEF.Sigma(tn_pick(ii),1,gmm_pick(ii)).^2).*COEF.SigmaSplit(tn_pick(ii),15)); % Natural logarithmic base
        Sinter(:,ii) = sqrt((COEF.Sigma(tn_pick(ii),1,gmm_pick(ii)).^2).*COEF.SigmaSplit(tn_pick(ii),14)); % Natural logarithmic base
    end
        GMmed(:,:) = (10.^immedtmp(:,:))/981; % (g)
       
        
end 


  



