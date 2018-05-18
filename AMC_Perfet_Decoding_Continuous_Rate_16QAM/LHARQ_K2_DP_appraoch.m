clc;
clear;
%clf;


SNRdB=[-5:30];
% correlation coefficient
rho=0.8 %(besselj(0,2*pi*0.2))

% Legend integral 


erg=@(x,av_snr) log2(1+x).*exp(-x./av_snr)./av_snr;

tol_step=1e-2;  
tol_abs=1e-2;
lambdaMin=1e-1;    
lambdaMax=14;
Nb_itt=2e2;
K=2;

%% Computation of throughput
for ind_snr=15:length(SNRdB)
    ind_snr
    snr=10.^(SNRdB(ind_snr)./10);
    ergodic_capacity(ind_snr)=integral(@(x) erg(x,snr),0,inf);
    max_rate=1.5.*ergodic_capacity(ind_snr);
    Dx=snr.*30;
    
    [bon_lambda]=Brent(@(lambda) find_lambda_DP(lambda,snr,max_rate,K,rho),lambdaMin,lambdaMax,tol_step, tol_abs,Nb_itt)
    
    
   

    
end