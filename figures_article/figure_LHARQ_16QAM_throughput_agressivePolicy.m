clc;
clear;
%clf;

figure(1)

load Er_MI_16QAM.mat
colour_peach = [251 111 66] ./ 255;
colour_darkblue = [1 17 181] ./ 255;

plot(SNR,Er_MI_16QAM_fading,'r','linewidth',1.3)

hold on 
plot(10,10,'-sm','markerfacecolor','m','MarkerSize',7.5)
plot(10,10,'-^b','markerfacecolor','b','MarkerSize',7.5)
plot(10,10,'-vk','markerfacecolor','k','MarkerSize',7.5)
plot(10,10,'-og','markerfacecolor','g','MarkerSize',7.5)

hold on


load LHARQ_AMC_policy_16QAM_Perfect_Decoding_continuous_R_dop_0.05_rho_0.95_V_leger.mat

plot(SNRdB,throughput_harq(1,:),'-m','linewidth',1.2)
amrk1=throughput_harq(1,(3:5:end-2));

plot(SNRdB,throughput_harq(4,:),'-b','linewidth',1.2)
amrk2=throughput_harq(4,(5:5:end-2));

load LHARQ_translated_AMC_policy_16QAM_Perfect_Decoding_continuous_R_dop_0.05rho_0.95_V_leger.mat
plot(SNRdB,throughput_harq(4,:),'-k','linewidth',1.2)
amrk3=throughput_harq(4,(6:5:end-2));

load HARQ_Dropping_translated_AMC_policy_16QAM_K_4_Perfect_Decoding_continuous_R_dop_0.05_rho_0.95_V_LHARQ_Correct.mat
plot(SNRdB,throughput_lharq,'-g','linewidth',1.2)
amrk4=throughput_lharq((4:5:end-2));


plot(SNRdB(3:5:end-2),amrk1,'sm','markerfacecolor','m','MarkerSize',7.5)

plot(SNRdB(5:5:end-2),amrk2,'^b','markerfacecolor','b','MarkerSize',7.5)

plot(SNRdB(6:5:end-2),amrk3,'vk','markerfacecolor','k','MarkerSize',7.5)

plot(SNRdB(4:5:end-2),amrk4,'og','markerfacecolor','g','MarkerSize',7.5)



xlabel('SNR')
ylabel('throughput')

grid

legend('Erg','AMC','HARQ','HARQT','DroppingX','Location','SouthEast')
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




name_figure = 'Agressive_dropping_HARQ';
print('-depsc','-r300',name_figure)
