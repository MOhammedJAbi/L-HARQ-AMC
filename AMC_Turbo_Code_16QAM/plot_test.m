clc;
clear;
clf;
figure(5)

load Er_MI_16QAM.mat

hold on

plot(SNR,Er_MI_16QAM_fading,'r','linewidth',1.1)
plot(10,10,'-m','markerfacecolor','m','MarkerSize',7)
plot(10,10,'-b','markerfacecolor','b','MarkerSize',7)
plot(10,10,'-k','linewidth',1.2)
plot(10,10,'--k','linewidth',1.2)


load LHARQ_translated_AMC_policy_16QAM_turbo_code_K_1_dop_0.1_rho_0.8.mat
plot(SNRdB,throughput_amc,'m','linewidth',1.1)

 load HARQ_AMC_policy_16QAM_Turbo_code_dop_0.1rho_0.8167.mat
plot(SNRdB,AMC_eta(4,:),'-b','linewidth',1.1)

load LHARQ_translated_AMC_policy_16QAM_turbo_code_K_1_dop_0.05_rho_0.95.mat
plot(SNRdB,throughput_amc,'--m','linewidth',1.1)
 
load HARQ_AMC_policy_16QAM_turbo_code_dop_0.05_rho_0.95_K1234_opt.mat
plot(SNRdB,AMC_eta(4,:),'--b','linewidth',1.1)

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




name_figure = 'TC_HARQ';
print('-depsc','-r300',name_figure)








figure(2)

load Er_MI_16QAM.mat
plot(SNR,Er_MI_16QAM_fading,'k','linewidth',1.1)

hold on

 load LHARQ_translated_AMC_policy_16QAM_turbo_code_K_4_dop_0.05_rho_0.95_epsilon_0.1_smallSet.mat
 plot(SNRdB,throughput_lharq,'r','linewidth',1.1)

 load LHARQ_translated_AMC_policy_16QAM_turbo_code_K_2_dop_0.05_rho_0.95_epsilon_0.1_smallSet.mat
 plot(SNRdB,throughput_lharq,'-.g','linewidth',1.1)
 
 
load LHARQ_translated_AMC_policy_16QAM_turbo_code_K_1_dop_0.05_rho_0.95.mat
plot(SNRdB,throughput_amc,'-m','linewidth',1.1)
 
 load HARQ_AMC_policy_16QAM_turbo_code_dop_0.05_rho_0.95_K1234_opt.mat
plot(SNRdB,AMC_eta(4,:),'-.b','linewidth',1.1)
%legend('Ergodic','LHARQ,K=4,Optimal Delta','LHARQ,K=2,Optimal Delta','AMC','HARQ,K=4,Delta=1')
grid
legend(' Erg',' LH4O',' LH2O',' AM',' HA','Location','SouthEast')
axis([-5 30 0 4])

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




name_figure = 'TB_LHARQ_rho95';
print('-depsc','-r300',name_figure)


figure(3)

load Er_MI_16QAM.mat
plot(SNR,Er_MI_16QAM_fading,'k','linewidth',1.1)

hold on

 load LHARQ_translated_AMC_policy_16QAM_turbo_code_K_4_dop_0.05_rho_0.95_epsilon_0.1_smallSet.mat
 plot(SNRdB,throughput_lharq,'r','linewidth',1.1)
 
 plot(SNRdB, throughput_lharq_delta(1,:),'-.k','linewidth',1.1)


 load LHARQ_translated_AMC_policy_16QAM_turbo_code_K_2_dop_0.05_rho_0.95_epsilon_0.1_smallSet.mat
 plot(SNRdB,throughput_lharq,'-.g','linewidth',1.1)
 plot(SNRdB, throughput_lharq_delta(1,:),'-.k','linewidth',1.1)
 
 
load LHARQ_translated_AMC_policy_16QAM_turbo_code_K_1_dop_0.05_rho_0.95.mat
plot(SNRdB,throughput_amc,'-m','linewidth',1.1)
 
 load HARQ_AMC_policy_16QAM_turbo_code_dop_0.05_rho_0.95_K1234_opt.mat
plot(SNRdB,AMC_eta(4,:),'-.b','linewidth',1.1)
legend('Ergodic','LHARQ,K=4,Optimal Delta','LHARQ,K=4,Delta=1','LHARQ,K=2,Optimal Delta','LHARQ,K=4,Delta=1','AMC','HARQ,K=4,Delta=1')

axis([-5 30 0 3.75])
grid

title('rho=0.95')


figure(1)

load Er_MI_16QAM.mat
plot(SNR,Er_MI_16QAM_fading,'k','linewidth',1.1)

hold on

 load LHARQ_translated_AMC_policy_16QAM_turbo_code_K_4_dop_0.1_rho_0.8_epsilon_0.1_smallSet.mat
 plot(SNRdB,throughput_lharq,'r','linewidth',1.1)
% 
 load LHARQ_translated_AMC_policy_16QAM_turbo_code_K_2_dop_0.1_rho_0.8_epsilon_0.1_smallSet.mat
 plot(SNRdB,throughput_lharq,'-.g','linewidth',1.1)


load LHARQ_translated_AMC_policy_16QAM_turbo_code_K_1_dop_0.1_rho_0.8.mat
plot(SNRdB,throughput_amc,'m','linewidth',1.1)

 load HARQ_AMC_policy_16QAM_Turbo_code_dop_0.1rho_0.8167.mat
plot(SNRdB,AMC_eta(4,:),'-.b','linewidth',1.1)

grid
legend(' Erg',' LH4O',' LH2O',' AM',' HA','Location','SouthEast')
axis([-5 30 0 4])

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




name_figure = 'TB_LHARQ_rho80';
print('-depsc','-r300',name_figure)
