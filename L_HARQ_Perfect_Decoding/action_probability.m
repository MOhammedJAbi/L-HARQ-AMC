clear
clc

load TS_fixed_rate_var_snr_R_4_K_3.mat
%load fixed_rate_var_snr_R_V8.mat

% Pr1=zeros(1,length(SNR));
% Pr2=zeros(1,length(SNR));
% Pr3=zeros(1,length(SNR));
% Pr4=zeros(1,length(SNR));

% Pr1_montec=zeros(1,length(SNR));
% Pr2_montec=zeros(1,length(SNR));
% Pr3_montec=zeros(1,length(SNR));
% Pr4_montec=zeros(1,length(SNR));

for j=1:length(SNR)

    j

    [P_sta,P]=states_steady_probabilities(Policy_corre(j,:),M,R,SNR(j));  %(R(k),SNR(j),step,policy_sc_ad2(:,j,k));
    
    steady_state(j,:)=P_sta;
    
    
%     Pn=P_sta(1:(M+1));
%     Pin=sum(Pn);
%     P_ack(j)=P_sta(M+1)/Pin;
%     POLY=Policy_corre(j,1:(M+1));
%     P_drop(j)=(sum(Pn(find(POLY==0)))-P_sta(M+1))/Pin;
%     P_ret(j)=sum(Pn(find(POLY==2)))/Pin;
%     P_ts(j)=sum(Pn(find(POLY~=0&POLY~=2)))/Pin;
    
%     PP=sum(P_sta(1:2*(step+1)))+P(2*(step+1)+1,2*(step+1)+1)*P_sta(2*(step+1)+1)+P(2*(step+1)+2,2*(step+1)+1)*P_sta(2*(step+1)+2);
%     Pn=P_sta(1:step)+P_sta(step+2:2*step+1);
%     POLY=policy_sc_ad2(1:step,j,k);
%     Pr1(j)=1/PP*sum(Pn(find(POLY==0)));
%     Pr2(j)=1/PP*sum(Pn(find(POLY==2)));
%     Pr4(j)=1/PP*sum(Pn(find(POLY~=0&POLY~=2)));
%     Pr3(j)=1-Pr1(j)-Pr2(j)-Pr4(j);
%     PLOLI=policy_sc_ad2(1:step+1,j,k);
%     [T,ac_ret,ac_sc,ac_drop,ac_no_need]=MDP_evaluation_montecarlo(R(k),PLOLI,SNR(j));
%     Pr1_montec(j)=ac_drop;
%     Pr2_montec(j)=ac_ret;
%     Pr3_montec(j)=ac_no_need;
%     Pr4_montec(j)=ac_sc;
end


save(['TS_fixed_rate_var_snr_steady_state_R_', num2str(R),'_K_', num2str(K),'.mat'])


