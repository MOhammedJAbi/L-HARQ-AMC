clear;
clc;

dop=0.1
rho=(besselj(0,2*pi*dop)).^2;
load ['LHARQ_AMC_policy_Perfect_Decoding_continuous_R_dop_',num2str(dop),'_rho_',num2str(rho),'.mat']

save ['LHARQ_AMC_policy_Perfect_Decoding_continuous_R_dop_',num2str(dop),'_rho_',num2str(rho),'.V_leger.mat'] dop rho SNRdB ...
    throughput_lharq throughput_harq outage snr_range policy

