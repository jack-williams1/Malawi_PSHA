%                                                                  %
% Maximum likelihood estimation of magnitude-recurrence parameters %
%                                                                  %

function GRpara = GRrelation_MLEWeichert_EQMATca(count,mag,dm,time)

ini_beta = 1.5;

Sum_countmag = sum(count.*mag);
N_count = sum(count);

% Iteration process starts
Error = ini_beta;
while Error >= 0.000001;
    
    Sum_exp = sum(exp(-ini_beta*mag));
    Sum_texp = sum(time.*exp(-ini_beta*mag));
    Sum_tmagexp = sum(time.*mag.*exp(-ini_beta*mag));
    Sum_tmag2exp = sum(time.*(mag.^2).*exp(-ini_beta*mag));
    
    d2LdB2 = N_count*((Sum_tmagexp/Sum_texp)^2-Sum_tmag2exp/Sum_texp);
    dLdB = Sum_tmagexp/Sum_texp*N_count - Sum_countmag;
    
    Mean_Beta = (ini_beta -dLdB/d2LdB2);
    Std_Beta = sqrt(-1/d2LdB2); 
    Std_Beta2 = sqrt(abs(Sum_texp^2/(N_count*(Sum_tmagexp^2-Sum_texp*Sum_tmag2exp))));
    
    Mean_B = Mean_Beta/log(10);
    Std_B = Std_Beta/log(10);
    Std_B2 = Std_Beta2/log(10);
    
    Mean_Na = N_count*Sum_exp/Sum_texp;
    Std_Na = sqrt(Mean_Na/N_count);
    
    F_N5 = N_count*Sum_exp/Sum_texp*exp(-Mean_Beta*(5-(min(mag)-dm/2)));
    F_N0 = N_count*Sum_exp/Sum_texp*exp(Mean_Beta*(min(mag)-dm/2));
    Std_F_N5 = F_N5/sqrt(N_count);
    Std_F_N0 = F_N0/sqrt(N_count);
    
    Error = abs(Mean_Beta-ini_beta);
    ini_beta = Mean_Beta;
    
end

GRpara = [Mean_B Std_B F_N0 Std_F_N0];


