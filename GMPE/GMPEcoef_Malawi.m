%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%   Coefficients for ground motion prediction equations for Malawi   %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

function COEF = GMPEcoef_Malawi(GMPEindex,Tn)

num_Tn = length(Tn);% Number of vibration periods assessed

if GMPEindex == 1
    
    load GMPEcoef_Malawi BA08
    
    % Mw: 5.0 - 8 and Rjb < 200 km
    % PGA and PSA (g) and PGV (cm/s) (cm) - geometric mean
    % ln(Y) = FM + FD + FS %So ground shaking a function of focal mechanism
    % (FM), Focal Distance (FD), and site geology (FS)
    % FM = e1*U + e2*SS + e3*NS + e4*RS + e5*min(M-Mh,0) + e6*(min(M-Mh,0))^2 + e7*max(M-Mh,0)
    % FD = (c1 + c2*(M-Mref))*ln(sqrt(Rjb^2+h^2)/rref) + c3*(sqrt(Rjb^2+h^2)-rref)
    % FS = blin*ln(max(min(Vs30,1300),180)/760) + bnl*ln(0.06/0.1) + c*(ln(max(pga4nl,0.03)/0.03))^2 + d*(ln(max(pga4nl,0.03)/0.03))^3 for pga4nl <= 0.09
    % FS = blin*ln(max(min(Vs30,1300),180)/760) + bnl*ln(pga4nl/0.1) for pga4nl > 0.09
    % bnl = (b1-b2)*ln(min(max(Vs30,180),300)/300)/ln(180/300) + b2*ln(max(min(Vs30,760),300)/760)/ln(300/760) 
    % U: unspecified mechanism; SS: strike-slip mechanism; NS: normal-slip mechanism; RS: reverse-slip mechanism
    % Coef = [e1 e2 e3 e4 e5 e6 e7 Mh c1 c2 c3 h blin b1 b2 sigma tau-U sigma-TU tau-M sigma-TM]
    % Maximum T is 10.0 s        
    
    COEF = zeros(num_Tn,20);
       
    for i = 1:num_Tn
        
        if Tn(i) == 0 % PGA
            COEF(i,1:20) = BA08(2,2:21);
        elseif Tn(i) > 0 % PSA
            for j = 1:20
                COEF(i,j) = interp1q(log(BA08(3:23,1)),BA08(3:23,1+j),log(Tn(i))); 
            end
        end
        
    end    
    
elseif GMPEindex == 2
    
    load GMPEcoef_Malawi CB08
    
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
    
    COEF = zeros(num_Tn,22);
       
    for i = 1:num_Tn
        
        if Tn(i) == 0 % PGA
            COEF(i,1:22) = CB08(3,2:23);
        elseif Tn(i) > 0 % PSA
            for j = 1:22
                COEF(i,j) = interp1q(log(CB08(4:24,1)),CB08(4:24,1+j),log(Tn(i))); 
            end
        end
        
    end
    
elseif GMPEindex == 3
    
    load GMPEcoef_Malawi CY08
    
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
    
    COEF = zeros(num_Tn,28);
       
    for i = 1:num_Tn
        
        if Tn(i) == 0 % PGA
            COEF(i,1:28) = CY08(2,2:29);
        elseif Tn(i) > 0 % PSA
            for j = 1:28
                COEF(i,j) = interp1q(log(CY08(3:24,1)),CY08(3:24,1+j),log(Tn(i))); 
            end
        end
        
    end
        
elseif GMPEindex == 4
    
    load GMPEcoef_Malawi AB10
    
    % log10(PGV/PGA/PSA) = b1 + b2*M + b3*M^2 + (b4+b5*M)*log10(sqrt(Rjb^2+b6^2)) + b7*Ss + b8*Sa + b9*Fn + b10*Fr
    % Ss = 1: soft soil site index (Vs30 < 360 m/s)
    % Sa = 1: firm soil site index (360 < Vs30 < 750 m/s)
    % Fn = 1: normal faulting index
    % Fr = 1: reverse faulting index
    % Constant sigma model
    % Coef = [b1 b2 b3 b4 b5 b6 b7 b8 b9 b10 sigma-intra sigma-inter sigma-total]
    % Maximum T is 3.0 s    
    
    COEF = zeros(num_Tn,13);
   
    for i = 1:num_Tn
        
        if Tn(i) == 0 % PGA
            COEF(i,1:13) = AB10(2,2:14);
        elseif Tn(i) > 0 % PSA
            for j = 1:13
                COEF(i,j) = interp1q(log(AB10(3:62,1)),AB10(3:62,1+j),log(Tn(i))); 
            end
        end
        
    end
    
elseif GMPEindex == 5
   
    load GMPEcoef_Malawi BSSA14
    
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

    COEF = zeros(num_Tn,31);
       
    for i = 1:num_Tn
        
        if Tn(i) == 0 % PGA
            COEF(i,1:31) = BSSA14(2,2:32);
        elseif Tn(i) > 0 % PSA
            for j = 1:31
                COEF(i,j) = interp1q(log(BSSA14(3:107,1)),BSSA14(3:107,1+j),log(Tn(i))); 
            end
        end
        
    end    
    
elseif GMPEindex == 6 || GMPEindex == 7 || GMPEindex == 8
    
    if GMPEindex == 6
        load GMPEcoef_Malawi ASB14_Rjb
        ASB14 = ASB14_Rjb;
    elseif GMPEindex == 7
        load GMPEcoef_Malawi ASB14_Repi
        ASB14 = ASB14_Repi;
    elseif GMPEindex == 8
        load GMPEcoef_Malawi ASB14_Rhypo
        ASB14 = ASB14_Rhypo;
    end
    
    % ln(Y) = ln(Yref(Mw,R,Fn,Fr)) + ln(S(Vs,PGAref)) + epsilon*sigma	
    % ln(Yref) = a1 + 0.0029*(Mw-6.75) + a3*(8.5-Mw)^2 + (a4+0.2529*(Mw-6.75))*ln(sqrt(R^2+7.5^2)) + a8*Fn + a9*Fr   for Mw <= 6.75	
    % ln(Yref) = a1 - 0.5096*(Mw-6.75) + a3*(8.5-Mw)^2 + (a4+0.2529*(Mw-6.75))*ln(sqrt(R^2+7.5^2)) + a8*Fn + a9*Fr   for Mw > 6.75	
    % ln(S) = b1*ln(Vs/750) + b2*ln([PGAref+2.5*(Vs/750)^3.2]/[(PGAref+2.5)*(Vs/750)^3.2])   for Vs <= 750	
    % ln(S) = b1*ln(min(Vs,1000)/750)   for Vs > 750	
    % PGA and PSA (g) and PGV (cm/s)	
    % R can be Rjb, Repi, Rhypo	
    % Y is geometric mean	
    % Mw: 4.0 to 7.6 and Rjb: 0 to 200 km	
    % Coef = [a1 a3 a4 a8 a9 b1 b2 sigma1 sigma2 sigma]
    % T: 0 to 4.0 s	
    
    COEF = zeros(num_Tn,10);
       
    for i = 1:num_Tn
        
        if Tn(i) == 0 % PGA
            COEF(i,1:10) = ASB14(2,2:11);
        elseif Tn(i) > 0 % PSA
            for j = 1:10
                COEF(i,j) = interp1q(log(ASB14(3:64,1)),ASB14(3:64,1+j),log(Tn(i))); 
            end
        end
        
    end    
    
elseif GMPEindex == 9
    
    load GMPEcoef_Malawi CFVB15
    
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
    
    COEF = zeros(num_Tn,20);
       
    for i = 1:num_Tn
        
        if Tn(i) == 0 % PGA
            COEF(i,1:20) = CFVB15(2,2:21);
        elseif Tn(i) > 0 % PSA
            for j = 1:20
                COEF(i,j) = interp1q(log(CFVB15(3:210,1)),CFVB15(3:210,1+j),log(Tn(i))); 
            end
        end
        
    end    

elseif GMPEindex == 10
    
    load GMPEcoef_Malawi CY14
    
    %  1   ) Period 
    %  2- 6) c1, c1_a, c1_b, c1_c, c1_d 
    %  7- 8) c7, c7_b 
    %  9-10) c11, c11_b
    % 11   ) c2 
    % 12   ) c3 
    % 13   ) cn 
    % 14   ) cm 
    % 15   ) c4 
    % 16   ) c5 
    % 17   ) c6 
    % 18   ) c_HM 
    % 19   ) c4a 
    % 20   ) c_RB 
    % 21-23) c_g1, c_g2, c_g3 
    % 24-26) c8, c8_a, c8_b 
    % 27-29) c9, c9_a, c9_b 
    % 30-35) phi1, phi2, phi3, phi4, phi5, phi6 
    % 36-38) sigma1, sigma2, sigma3 
    % 39-40) tau1, tau2
    
    COEF = zeros(num_Tn,size(CY14,2));
    
    for i = 1:num_Tn
        
        [~,indx] = min(abs(Tn(i)-CY14(:,1))); % index coef for correct spectral acceleration Tn in CY14
        for j = 1:size(CY14,2)
            COEF(i,j) = CY14(indx,j);
        end
        
    end
    
elseif GMPEindex == 11 
    
    % Ground motion values are tabulated for various combinations of magnitude, hypocentral distance, vibration period, Vs30, and logic tree branches [see Kolaj et al. 2020].
    
    load GMPEcoef_Malawi CanadaSHM6_2015
    
    COEF = CanadaSHM6_2015;
    
elseif GMPEindex == 12 
    
    % Ground motion values are tabulated for various combinations of magnitude, rupture distance, vibration period, Vs30, and logic tree branches [see Kolaj et al. 2020].

    load GMPEcoef_Malawi CanadaSHM6_2020
    
    COEF = CanadaSHM6_2020;
    
end 
