clc;
clear;
%clf;

figure(1)

load Er_MI_16QAM.mat

plot(SNR,Er_MI_16QAM_fading,'r','linewidth',1.3)

hold on 
plot(10,10,'-sm','markerfacecolor','m','MarkerSize',7)
plot(10,10,'-^b','markerfacecolor','b','MarkerSize',7)
plot(10,10,'-k','linewidth',1.2)
plot(10,10,'--k','linewidth',1.2)
hold on

load LHARQ_AMC_policy_16QAM_Perfect_Decoding_continuous_R_dop_0.1_rho_0.8_V_leger.mat

plot(SNRdB,throughput_harq(1,:),'-m','linewidth',1.2)
amrk1=throughput_harq(1,(3:5:end-2));

plot(SNRdB,throughput_harq(4,:),'-b','linewidth',1.2)
amrk2=throughput_harq(4,(5:5:end-2));


plot(SNRdB(3:5:end-2),amrk1,'sm','markerfacecolor','m','MarkerSize',7)

plot(SNRdB(5:5:end-2),amrk2,'^b','markerfacecolor','b','MarkerSize',7)

%%

load LHARQ_AMC_policy_16QAM_Perfect_Decoding_continuous_R_dop_0.05_rho_0.95_V_leger.mat

plot(SNRdB,throughput_harq(1,:),'--m','linewidth',1.2)
amrk1=throughput_harq(1,(3:5:end-2));

plot(SNRdB,throughput_harq(4,:),'--b','linewidth',1.2)
amrk2=throughput_harq(4,(5:5:end-2));


plot(SNRdB(3:5:end-2),amrk1,'sm','markerfacecolor','m','MarkerSize',7)

plot(SNRdB(5:5:end-2),amrk2,'^b','markerfacecolor','b','MarkerSize',7)





xlabel('SNR')
ylabel('throughput')

grid

legend('Erg','AMC','HARQ','t=0.1XXXX','t=0.05','Location','SouthEast')
axis([-5 30 0 4])

set(gca,'fontsize',11)
ll=legend;
set(ll,'FontSize',15);
%set(gcf,'PaperPositionMode','manual');

set(gcf,'PaperUnits','centimeters');
lar=512;
lon=2*lar/(1.1+sqrt(5));
set(gcf,'Position',[384 874 lar lon]);
set(gcf,'PaperPosition',[1.19 16 1.19 16]);
set(gca,'fontName','Times')
set(gcf,'PaperPositionMode','auto');




name_figure = 'HARQ';
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