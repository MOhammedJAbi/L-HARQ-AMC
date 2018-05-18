clc;
clear;
clf;

figure(1)

% corrected L-HARQ, problem of renewal-reward K=2
load LHARQ_translated_AMC_policy_16QAM_K_2_Perfect_Decoding_continuous_R_dop_0.05_rho_0.95_V_LHARQ_Correct_Part1.mat


%L-HARQ throughput with optimal delta
optimal_delta_K_2=optimal_delta(1,(1:18));

load LHARQ_translated_AMC_policy_16QAM_K_2_Perfect_Decoding_continuous_R_dop_0.05_rho_0.95_V_LHARQ_Correct_Part2.mat

optimal_delta_K_2=[optimal_delta_K_2 optimal_delta(1,(19:end))];
%

% corrected L-HARQ, problem of renewal-reward K=4
load LHARQ_translated_AMC_policy_16QAM_K_4_Perfect_Decoding_continuous_R_dop_0.05_rho_0.95_V_LHARQ_Correct_Part1.mat

%L-HARQ throughput with AMC_policy
optimal_delta_K_4=optimal_delta(1,(1:18));

load LHARQ_translated_AMC_policy_16QAM_K_4_Perfect_Decoding_continuous_R_dop_0.05_rho_0.95_V_LHARQ_Correct_Part2.mat

optimal_delta_K_4=[optimal_delta_K_4 optimal_delta(1,(19:end))];
%

plot(SNRdB, 10.*log10(optimal_delta_K_2),'r',SNRdB,10.*log10(optimal_delta_K_4),'b')