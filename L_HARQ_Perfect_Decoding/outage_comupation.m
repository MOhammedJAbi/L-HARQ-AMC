clear;
clc;

load TS_fixed_rate_var_snr_R_4_K_2_MI_31_Alpha_16.mat

for j=1:length(SNR)
    j
    s=(M+1)*(M+3); %# of states
    
    P_k_1=Matrix_transition(Policy_corre(j,:),M,R,SNR(j),ind_hp);
    I=eye(size(P_k_1));
    K=inv(I-P_k_1);
    
    rew_of_re_k=zeros(s,1);
    
    for i=1:s
        rew_of_re_k(i)=Np_expectation(i,Policy_corre(j,i),R,M,SNR(j));
    end
    
    O=K*rew_of_re_k;
    zeta_p=O(ind_hp)/(sum(K(ind_hp,:)));
    Pout(j)=1-T_corre(j)/(R*zeta_p);
end

save(['TS_Outage_fixed_rate_var_snr_R_', num2str(R),'_K_', num2str(K),'_MI_', num2str(M),'_Alpha_', num2str(L),'.mat'])


