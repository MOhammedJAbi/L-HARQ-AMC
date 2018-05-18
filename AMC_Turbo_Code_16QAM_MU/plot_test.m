clc;
clear;
clf;

figure(1)
load MU_LHARQ_translated_AMC_policy_16QAM_turbo_code_K_1_dop_0.05_rho_0.95.mat

plot(L_set,throughput_amc,'r')

hold on

load MU_LHARQ_translated_AMC_policy_16QAM_turbo_code_K_2_dop_0.05_rho_0.95.mat

plot(L_set,throughput_lharq,'b')

legend('AMC','LHARQ,K=2,Delta=1')

title('tho=0.95')

figure(2)
load MU_LHARQ_translated_AMC_policy_16QAM_turbo_code_K_1_dop_0.1_rho_0.8.mat

plot(L_set,throughput_amc,'r')

hold on

load MU_LHARQ_translated_AMC_policy_16QAM_turbo_code_K_2_dop_0.1_rho_0.8.mat

plot(L_set,throughput_lharq,'b')

legend('AMC','LHARQ,K=2,Delta=1')

title('tho=0.8')



