clc;
clear;
%clf;

figure(1)

load Er_MI_16QAM.mat


plot(SNR,Er_MI_16QAM_fading,'k','linewidth',1.1)


hold on

load AMC_policy_16QAM_Perfect_Decoding_Discrete_R_dop_0.05rho_0.95156.mat
plot(SNRdB,throughput,'k','linewidth',1.1)


load LHARQ_translated_AMC_policy_16QAM_K_2_Perfect_Decoding_Discrete_R_dop_0.05_rho_0.95_V_LHARQ_Correct_Vprovisoir.mat

plot(SNRdB,throughput_lharq,'r')

hold on

load LHARQ_translated_AMC_policy_16QAM_K_4_Perfect_Decoding_Discrete_R_dop_0.05_rho_0.95_V_LHARQ_Correct_Vprovisoir.mat

plot(SNRdB,throughput_lharq,'-b')

load LHARQ_translated_AMC_policy_16QAM_K_2_Perfect_Decoding_Discrete_R_dop_0.05_rho_0.95_V_LHARQ_Correct_Vprovisoir_limit_D2.mat

plot(SNRdB,throughput_lharq,'-.k')

hold on

load LHARQ_translated_AMC_policy_16QAM_K_4_Perfect_Decoding_Discrete_R_dop_0.05_rho_0.95_V_LHARQ_Correct_Vprovisoir_limit_D2.mat

plot(SNRdB,throughput_lharq,'-.m')

% load LHARQ_translated_AMC_policy_16QAM_K_2_Perfect_Decoding_Discrete_R_dop_0.05_rho_0.95_V_LHARQ_Correct.mat
% 
% plot(SNRdB,throughput_lharq,'-.g')


axis([-5 30 0 4])


legend('Erg','LHARQ,K=2,large set D','LHARQ,K=4,large set D','LHARQ,K=2,small set D','LHARQ,K=4,small set D','Location','SouthEast')





% name_figure = 'L_HARQ';
% print('-depsc','-r300',name_figure)

% figure(2)
% 
% load LHARQ_outage_policy_16QAM_Perfect_Decoding_continuous_R_dop_0.05_rho_0.95_V_leger_V4.mat
% 
% semilogy(SNRdB,optimal_epsilon_lharq(1,:),'r',SNRdB,optimal_epsilon_lharq(2,:),'b',SNRdB,optimal_epsilon_lharq(3,:),'g',SNRdB,optimal_epsilon_lharq(4,:),'m')
% hold on
% semilogy(SNRdB,optimal_epsilon_lharq_p(1,:),'-.r',SNRdB,optimal_epsilon_lharq_p(2,:),'-.b',SNRdB,optimal_epsilon_lharq_p(3,:),'-.g',SNRdB,optimal_epsilon_lharq_p(4,:),'-.m')
% grid;
% figure(3)
% 
% ii=20;
% semilogx(snr_range(ii,:),squeeze(policy_lharq(1,ii,:)),'r',snr_range(ii,:),squeeze(policy_lharq(2,ii,:)),'b',snr_range(ii,:),squeeze(policy_lharq(3,ii,:)),'g',snr_range(ii,:),squeeze(policy_lharq(4,ii,:)),'m')
% %axis([0 100 0 4]);
% grid;