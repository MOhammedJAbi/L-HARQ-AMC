clc;
clear;
%clf;

figure(1)

load Er_MI_16QAM.mat


plot(SNR,Er_MI_16QAM_fading,'k','linewidth',1.1)


hold on

load AMC_policy_16QAM_Perfect_Decoding_Discrete_R_dop_0.05rho_0.95156.mat
plot(SNRdB,throughput,'k','linewidth',1.1)

load HARQ_AMC_policy_16QAM_K_2_Perfect_Decoding_Discrete_R_dop_0.05_rho_0.95.mat

plot(SNRdB,throughput_harq,'m')



load LHARQ_translated_AMC_policy_16QAM_K_2_Perfect_Decoding_Discrete_R_dop_0.05_rho_0.95_V_LHARQ_Correct.mat

plot(SNRdB,throughput_lharq,'-.r')

plot(SNRdB,throughput_lharq_delta(1,:),'r')

hold on

load VLHARQ_translated_AMC_policy_16QAM_K_2_Perfect_Decoding_Discrete_R_dop_0.05_rho_0.95.mat

plot(SNRdB,throughput_vlharq_delta(1,:),'-.g')

plot(SNRdB,throughput_vlharq(:),'--k')

legend(' Ergodic',' AMC','IRHARQ',' L-HARQ',' LHARQ-OP',' VL-HARQ',' VLHARQ-OP','Location','SouthEast')

