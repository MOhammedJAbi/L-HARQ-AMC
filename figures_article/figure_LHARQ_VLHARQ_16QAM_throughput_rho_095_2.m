clc;
clear;
%clf;
colour_peach = [251 111 66] ./ 255;
colour_darkblue = [1 17 181] ./ 255;
figure(1)

load Er_MI_16QAM.mat

plot(SNR,Er_MI_16QAM_fading,'r','linewidth',1.3)

hold on 
plot(10,10,'-sm','markerfacecolor','m','MarkerSize',7)
plot(10,10,'-^b','markerfacecolor','b','MarkerSize',7)
plot(10,10,'--og','markerfacecolor','g','MarkerSize',7)
plot(10,10,'--^','Color',colour_peach,'MarkerEdgeColor',colour_peach,'markerfacecolor',colour_peach,'MarkerSize',7)
plot(10,10,'-.o','Color',colour_darkblue,'MarkerEdgeColor',colour_darkblue,'markerfacecolor',colour_darkblue,'MarkerSize',7)
plot(10,10,'-.^k','markerfacecolor','k','MarkerSize',7)
hold on

% corrected L-HARQ, problem of renewal-reward K=2
load LHARQ_translated_AMC_policy_16QAM_K_2_Perfect_Decoding_continuous_R_dop_0.05_rho_0.95_V_LHARQ_Correct_Part1.mat

%L-HARQ throughput with AMC_policy
throughput_lharq_delta_1_K_2=throughput_lharq_delta(1,(1:18));

load LHARQ_translated_AMC_policy_16QAM_K_2_Perfect_Decoding_continuous_R_dop_0.05_rho_0.95_V_LHARQ_Correct_Part2.mat
throughput_lharq_delta_1_K_2=[throughput_lharq_delta_1_K_2 throughput_lharq_delta(1,(19:end))];
%

% corrected L-HARQ, problem of renewal-reward K=4
load LHARQ_translated_AMC_policy_16QAM_K_4_Perfect_Decoding_continuous_R_dop_0.05_rho_0.95_V_LHARQ_Correct_Part1.mat
throughput_lharq_delta_1_K_4=throughput_lharq_delta(1,(1:18));
%L-HARQ throughput with optimal delta

load LHARQ_translated_AMC_policy_16QAM_K_4_Perfect_Decoding_continuous_R_dop_0.05_rho_0.95_V_LHARQ_Correct_Part2.mat
throughput_lharq_delta_1_K_4=[throughput_lharq_delta_1_K_4 throughput_lharq_delta(1,(19:end))];
%


load LHARQ_AMC_policy_16QAM_Perfect_Decoding_continuous_R_dop_0.05_rho_0.95_V_leger.mat
plot(SNRdB,throughput_harq(1,:),'-m','linewidth',1.2)
amrk1=throughput_harq(1,(3:5:end-2));

plot(SNRdB,throughput_harq(4,:),'-b','linewidth',1.2)
amrk2=throughput_harq(4,(4:5:end-2));


plot(SNRdB(3:5:end-2),amrk1,'sm','markerfacecolor','m','MarkerSize',7)
plot(SNRdB(4:5:end-2),amrk2,'^b','markerfacecolor','b','MarkerSize',7)

%%

plot(SNRdB,throughput_lharq_delta_1_K_2,'--g','linewidth',1.2)
amrk1=throughput_lharq_delta_1_K_2(3:5:end-2);
plot(SNRdB,throughput_lharq_delta_1_K_4,'--','Color',colour_peach,'linewidth',1.2)
amrk2=throughput_lharq_delta_1_K_4(4:5:end-2);


load VLHARQ_translated_AMC_policy_16QAM_K_2_Perfect_Decoding_continous_R_dop_0.05_rho_0.95.mat
plot(SNRdB,throughput_vlharq_delta(1,:),'-.','Color',colour_darkblue,'linewidth',1.2)
amrk11=throughput_vlharq_delta(1,6:5:end-2);

load VLHARQ_translated_AMC_policy_16QAM_K_4_Perfect_Decoding_continous_R_dop_0.05_rho_0.95.mat
plot(SNRdB,throughput_vlharq_delta(1,:),'-.k','linewidth',1.2)
amrk22=throughput_vlharq_delta(1,(7:5:end-2));

plot(SNRdB(3:5:end-2),amrk1,'og','markerfacecolor','g','MarkerSize',7)

plot(SNRdB(4:5:end-2),amrk2,'^','MarkerEdgeColor',colour_peach,'markerfacecolor',colour_peach,'MarkerSize',7)

plot(SNRdB(6:5:end-2),amrk11,'o','MarkerEdgeColor',colour_darkblue,'markerfacecolor',colour_darkblue,'MarkerSize',7)

plot(SNRdB(7:5:end-2),amrk22,'^k','markerfacecolor','k','MarkerSize',7)

% load Joint_K_4_R_max_8_Pas_0.1.mat
% plot(SNRdB,T_corre,'k')




xlabel('SNR')
ylabel('throughput')

grid

legend(' Erg',' AM',' HA',' LH2',' LH4',' VLH2',' VLH4','Location','SouthEast')
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




name_figure = 'VL-LHARQ_rho95';
print('-depsc','-r300',name_figure)

