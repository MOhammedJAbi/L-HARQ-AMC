% optimization of modulation in each avarage snr

clc;
clear;
%clf;
colour_peach = [251 111 66] ./ 255;
colour_darkblue = [1 17 181] ./ 255;
figure(1)

load Er_MI_64QAM.mat

plot(SNR,Er_MI_64QAM_fading,'r','linewidth',1.3)

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
%L-HARQ throughput with optimal delta
throughput_lharq_K_2=throughput_lharq(1,(1:18));

load LHARQ_translated_AMC_policy_16QAM_K_2_Perfect_Decoding_continuous_R_dop_0.05_rho_0.95_V_LHARQ_Correct_Part2.mat

throughput_lharq_delta_1_K_2=[throughput_lharq_delta_1_K_2 throughput_lharq_delta(1,(19:end))];
throughput_lharq_K_2=[throughput_lharq_K_2 throughput_lharq(1,(19:end))];
%

% corrected L-HARQ, problem of renewal-reward K=4
load LHARQ_translated_AMC_policy_16QAM_K_4_Perfect_Decoding_continuous_R_dop_0.05_rho_0.95_V_LHARQ_Correct_Part1.mat

%L-HARQ throughput with AMC_policy
throughput_lharq_delta_1_K_4=throughput_lharq_delta(1,(1:18));
%L-HARQ throughput with optimal delta
throughput_lharq_K_4=throughput_lharq(1,(1:18));

load LHARQ_translated_AMC_policy_16QAM_K_4_Perfect_Decoding_continuous_R_dop_0.05_rho_0.95_V_LHARQ_Correct_Part2.mat

throughput_lharq_delta_1_K_4=[throughput_lharq_delta_1_K_4 throughput_lharq_delta(1,(19:end))];
throughput_lharq_K_4=[throughput_lharq_K_4 throughput_lharq(1,(19:end))];
%


load LHARQ_AMC_policy_16QAM_Perfect_Decoding_continuous_R_dop_0.05_rho_0.95_V_leger.mat
throughput_amc_16=[throughput_harq(1,:), zeros(1,10)];
throughput_harq_16=[throughput_harq(4,:), zeros(1,10)];


load LHARQ_AMC_policy_64QAM_Perfect_Decoding_continuous_R_dop_0.05rho_0.95_V_leger.mat
throughput_amc=max(throughput_harq(1,:),throughput_amc_16);
plot(SNRdB,throughput_amc,'-m','linewidth',1.2)
amrk1=throughput_amc(1,(3:5:end-2));

throughput_harq_16_64=max(throughput_harq_16,throughput_harq(4,:));
plot(SNRdB,throughput_harq_16_64,'-b','linewidth',1.2)
amrk2=throughput_harq_16_64(1,(4:5:end-2));


plot(SNRdB(3:5:end-2),amrk1,'sm','markerfacecolor','m','MarkerSize',7)

plot(SNRdB(4:5:end-2),amrk2,'^b','markerfacecolor','b','MarkerSize',7)

%%

load LHARQ_translated_AMC_policy_64QAM_K_2_Perfect_Decoding_continuous_R_dop_0.05_rho_0.95_V_LHARQ_Correct.mat
throughput_lharq_delta_1_K_2=[throughput_lharq_delta_1_K_2 zeros(1,10)];
throughput_lharq_delta_1_K_2=max(throughput_lharq_delta_1_K_2,throughput_lharq_delta(1,:));
plot(SNRdB,throughput_lharq_delta_1_K_2,'--g','linewidth',1.2)
amrk1=throughput_lharq_delta_1_K_2(3:5:end-2);

load LHARQ_translated_AMC_policy_64QAM_K_4_Perfect_Decoding_continuous_R_dop_0.05_rho_0.95_V_LHARQ_Correct.mat
throughput_lharq_delta_1_K_4=[throughput_lharq_delta_1_K_4 zeros(1,10)];
throughput_lharq_delta_1_K_4=max(throughput_lharq_delta_1_K_4,throughput_lharq_delta(1,:));
plot(SNRdB,throughput_lharq_delta_1_K_4,'--','Color',colour_peach,'linewidth',1.2)
amrk2=throughput_lharq_delta_1_K_4(4:5:end-2);


load LHARQ_translated_AMC_policy_64QAM_K_2_Perfect_Decoding_continuous_R_dop_0.05_rho_0.95_V_LHARQ_Correct.mat
throughput_lharq_K_2=[throughput_lharq_K_2 zeros(1,10)];
throughput_lharq_K_2=max(throughput_lharq_K_2,throughput_lharq);
plot(SNRdB,throughput_lharq_K_2,'-.','Color',colour_darkblue,'linewidth',1.2)
amrk11=throughput_lharq_K_2(6:5:end-2);

load LHARQ_translated_AMC_policy_64QAM_K_4_Perfect_Decoding_continuous_R_dop_0.05_rho_0.95_V_LHARQ_Correct.mat
throughput_lharq_K_4=[throughput_lharq_K_4 zeros(1,10)];
throughput_lharq_K_4=max(throughput_lharq_K_4,throughput_lharq);
plot(SNRdB,throughput_lharq_K_4,'-.k','linewidth',1.2)
amrk22=throughput_lharq_K_4(7:5:end-2);

plot(SNRdB(3:5:end-2),amrk1,'og','markerfacecolor','g','MarkerSize',7)

plot(SNRdB(4:5:end-2),amrk2,'^','MarkerEdgeColor',colour_peach,'markerfacecolor',colour_peach,'MarkerSize',7)

plot(SNRdB(6:5:end-2),amrk11,'o','MarkerEdgeColor',colour_darkblue,'markerfacecolor',colour_darkblue,'MarkerSize',7)

plot(SNRdB(7:5:end-2),amrk22,'^k','markerfacecolor','k','MarkerSize',7)




xlabel('SNR')
ylabel('throughput')

grid

legend(' Erg',' AM',' HA',' LH2',' LH4',' LH2O',' LH4O','Location','SouthEast')
axis([-5 40 0 6])

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




name_figure = 'LHARQ_rho95_16_64_QAM';
print('-depsc','-r300',name_figure)

