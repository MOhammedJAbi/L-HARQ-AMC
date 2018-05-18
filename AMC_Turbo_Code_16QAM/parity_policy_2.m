function [parity_policy] = parity_policy_2(epsilon_Set)

epsilon_Set;
%% R=1.5
load WEP_Turbo_16QAM_rate1.5poincenage6_3gpp.mat
parity_policy(1,:)=Rate_parity_Set(2:end-1);

%% R=2.25
SNRLin=0;
WEP_set=0;
load WEP_Turbo_16QAM_rate2.25poincenage6_3gpp.mat
parity_policy(2,:)=Rate_parity_Set(2:end-1);


%% R=3
SNRLin=0;
WEP_set=0;
load WEP_Turbo_16QAM_rate3poincenage6_3gpp.mat
parity_policy(3,:)=Rate_parity_Set(2:end-1);



%% R=3.75
SNRLin=0;
WEP_set=0;
load WEP_Turbo_16QAM_rate3.75poincenage6_3gpp.mat
parity_policy(4,:)=Rate_parity_Set(2:end-1);

