clear;
clc;

load fixed_rate_var_snr_R_V12.mat

plot(SNR,T_corre,'g')

hold on

load fixed_rate_var_snr_R_V8.mat

plot(SNR,T_corre,'-.r')

hold on

load fixed_rate_var_snr_R_V2.mat

plot(SNR,T_corre,'-k')

hold on

load constant_rate.mat

plot(snrdB,throughput(3,:,1),'b',snrdB,throughput(6,:,1),'-.m')

%axis([0 30 0 5])

grid