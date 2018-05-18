clc;
clear;
clf;

figure(1)



load LHARQ_translated_AMC_policy_16_64QAM_K_2_Perfect_Decoding_continuous_R_dop_0.05_rho_0.95_V_LHARQ_Correct.mat

ind=26;
SNRdB(ind)

plot(10.*log10(snr_range(ind,:)),policy(ind,:),'-m','linewidth',1.2)


hold on
plot(10.*log10(snr_range(ind,:)./optimal_delta(ind)),policy(ind,:),'--b','linewidth',1.2)

load LHARQ_translated_AMC_policy_16_64QAM_K_4_Perfect_Decoding_continuous_R_dop_0.05_rho_0.95_V_LHARQ_Correct.mat

ind=26;

plot(10.*log10(snr_range(ind,:)./optimal_delta(ind)),policy(ind,:),'-.k','linewidth',1.2)


xlabel('SNR')
ylabel('throughput')

grid

legend(' AMCC',' LHAQR2',' LHARQ4XXX','Location','NorthWest')
axis([-5 25 1.5 6])


set(gca,'fontsize',11)
ll=legend;
set(ll,'FontSize',16);
%set(gcf,'PaperPositionMode','manual');

set(gcf,'PaperUnits','centimeters');
lar=512;
lon=2*lar/(1.1+sqrt(5));
set(gcf,'Position',[384 874 lar lon]);
set(gcf,'PaperPosition',[1.19 16 1.19 16]);
set(gca,'fontName','Times')
set(gcf,'PaperPositionMode','auto');




name_figure = 'policy_LHARQ_rho95';
print('-depsc','-r300',name_figure)

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