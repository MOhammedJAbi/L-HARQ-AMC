clc;
clear;
clf('reset')

figure(1)

load LHARQ_AMC_DP_policy_Perfect_Decoding_continuous_R_dop_0.1_rho_0.8_V3.mat
plot(SNRdB,throughput_lharq(2,:),'b')

hold on

load LHARQ_AMC_DP_policy_Perfect_Decoding_continuous_R_dop_0.1_rho_0.8_V2.mat
plot(SNRdB,throughput_lharq(2,:),':b')

load LHARQ_AMC_policy_Perfect_Decoding_continuous_R_dop_0.1_rho_0.8_V_leger.mat

plot(SNRdB,throughput_lharq(1,:),'-.r',SNRdB,throughput_lharq(2,:),'-.b')

grid;
legend('DP3-LHARQ,K=2','DP2-LHARQ,K=2','AMC','AMC-LHARQ,K=2')







figure(2)
load LHARQ_outage_policy_Perfect_Decoding_continuous_R_dop_0.1_rho_0.8_V_leger.mat

plot(SNRdB,throughput_lharq(1,:),'r',SNRdB,throughput_lharq(2,:),'b',SNRdB,throughput_lharq(3,:),'g',SNRdB,throughput_lharq(4,:),'m')

hold on

load LHARQ_AMC_policy_Perfect_Decoding_continuous_R_dop_0.075_rho_0.9_V_leger.mat

plot(SNRdB,throughput_lharq(1,:),'-.r',SNRdB,throughput_lharq(2,:),'-.b',SNRdB,throughput_lharq(3,:),'-.g',SNRdB,throughput_lharq(4,:),'-.m')

hold on

% load Er_MI_16QAM.mat

%plot(SNRdB,ergodic_capacity,'k')

legend('K=1','K=2','K=3','K=4')
%axis([-5 30 0 4]);
grid;

figure(3)

load LHARQ_outage_policy_Perfect_Decoding_continuous_R_dop_0.1_rho_0.8_V_leger_V2.mat

semilogy(SNRdB,optimal_epsilon_lharq(1,:),'r',SNRdB,optimal_epsilon_lharq(2,:),'b',SNRdB,optimal_epsilon_lharq(3,:),'g',SNRdB,optimal_epsilon_lharq(4,:),'m')
hold on
load LHARQ_AMC_policy_Perfect_Decoding_continuous_R_dop_0.1_rho_0.8_V_leger.mat
semilogy(SNRdB,outage,'k')

% hold on
%  semilogy(SNRdB,optimal_epsilon_lharq_p(1,:),'-.r',SNRdB,optimal_epsilon_lharq_p(2,:),'-.b',SNRdB,optimal_epsilon_lharq_p(3,:),'-.g',SNRdB,optimal_epsilon_lharq_p(4,:),'-.m')
% grid;
legend('K=1','K=2','K=3','K=4')
figure(4)

ii=2;
plot(snr_range(ii,:),squeeze(policy_lharq(1,ii,:)),'r',snr_range(ii,:),squeeze(policy_lharq(2,ii,:)),'b',snr_range(ii,:),squeeze(policy_lharq(3,ii,:)),'g',snr_range(ii,:),squeeze(policy_lharq(4,ii,:)),'m')
%axis([0 100 0 4]);
legend('K=1','K=2','K=3','K=4')
grid;

figure(5)
load LHARQ_AMC_DP_policy_Perfect_Decoding_continuous_R_dop_0.1_rho_0.8_V3.mat
plot(SNRdB,optimal_lambda_lharq,'b')

hold on

load LHARQ_AMC_DP_policy_Perfect_Decoding_continuous_R_dop_0.1_rho_0.8_V2.mat
plot(SNRdB,optimal_lambda_lharq,':b')

legend('DP3-LHARQ,K=2','DP2-LHARQ,K=2')