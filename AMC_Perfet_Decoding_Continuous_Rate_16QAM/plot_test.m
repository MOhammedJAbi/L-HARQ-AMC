clc;
clear;
%clf;

figure(1)
% load LHARQ_outage_policy_16QAM_Perfect_Decoding_continuous_R_dop_0.05_rho_0.95_V_leger_V4.mat
% 
% plot(SNRdB,throughput_lharq(1,:),'r',SNRdB,throughput_lharq(2,:),'b',SNRdB,throughput_lharq(3,:),'g',SNRdB,throughput_lharq(4,:),'m')
% 
% 
% hold on

load Er_MI_16QAM.mat


plot(SNR,Er_MI_16QAM_fading,'k','linewidth',1.1)


hold on

load LHARQ_AMC_policy_16QAM_Perfect_Decoding_continuous_R_dop_0.075_rho_0.9_V_leger.mat

plot(SNRdB,throughput_lharq(1,:),'-r','linewidth',1.1)
amrk1=throughput_lharq(1,(3:3:end));

plot(SNRdB,throughput_lharq(2,:),'-g','linewidth',1.1)


plot(SNRdB,throughput_lharq(4,:),'-b','linewidth',1.1)

load LHARQ_translated_AMC_policy_16QAM_Perfect_Decoding_continuous_R_dop_0.075rho_0.9_V_leger.mat

plot(SNRdB,throughput_lharq(1,:),'-.k','linewidth',1.1)
amrk1=throughput_lharq(1,(3:3:end));

plot(SNRdB,throughput_harq(2,:),'-.k','linewidth',1.1)


plot(SNRdB,throughput_harq(4,:),'-.k','linewidth',1.1)
axis([-5 30 0 4])


legend('Erg','K=1','K=2','K=4','Location','SouthEast')

figure(2)

semilogy(SNRdB,delta_harq(2,:),'r',SNRdB,delta_harq(4,:),'b')



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