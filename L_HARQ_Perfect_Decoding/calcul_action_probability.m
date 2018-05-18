clear;
clc;

load TS_fixed_rate_var_snr_steady_state_R_4_K_3.mat

for j=1:length(SNR)

    j
    
    P_sta=steady_state(j,:);
    %%
    %G1
    P_G1=P_sta(1:(M+1));
    P_G1_norm=sum(P_G1);
    P_G1_ack(j)=P_sta(M+1)/P_G1_norm;
    POLY_G1=Policy_corre(j,1:(M+1));
    P_G1_drop(j)=(sum(P_G1(find(POLY_G1==0)))-P_G1(M+1))/P_G1_norm;
    P_G1_ret(j)=sum(P_G1(find(POLY_G1==2)))/P_G1_norm;
    P_G1_ts(j)=sum(P_G1(find(POLY_G1~=0&POLY_G1~=2)))/P_G1_norm;
    %%
    %G2
    P_G2=P_sta((M+1)+1:(M+1)*(M+2));
    P_G2_norm=sum(P_G2);
    P_G2_ack(j)=P_G2(end)/P_G2_norm;
    POLY_G2=Policy_corre(j,(M+1)+1:(M+1)*(M+2));
    %->G2&G3
    POLY_G2_case1=Policy_corre(j,(M+1)*(M+1)+1:(M+1)*(M+2)-1);
    P_G2_G2(j)=sum(P_G2(find(POLY_G2_case1~=0&POLY_G2_case1~=2)+M*(M+1)))/P_G2_norm;
    P_G2_G3_ret(j)=sum(P_G2(find(POLY_G2_case1==2)+M*(M+1)))/P_G2_norm;
    %->G4&G5
    POLY_G2_case2=Policy_corre(j,(M+1)+1:(M+1)*(M+1));
    for k=1:M
        POLY_G2_case2(k*(M+1))=3;
    end
    P_G2_G4(j)=sum(P_G2(find(POLY_G2_case2~=0&POLY_G2_case2~=2&POLY_G2_case2~=3)))/P_G2_norm;
    P_G2_G5_ret(j)=sum(P_G2(find(POLY_G2_case2==2)))/P_G2_norm;
    P_G2_G3_drop(j)=sum(P_G2(find(POLY_G2_case2==0)))/P_G2_norm;
    %->G5&G6
    POLY_G2_case3=Policy_corre(j,(M+1)*2:(M+1):(M+1)*(M+1));
    P_G2_G5_ts(j)=sum(P_G2((find(POLY_G2_case3~=0&POLY_G2_case3~=2))*(M+1)))/P_G2_norm;
    P_G2_G6_ret(j)=sum(P_G2(find(POLY_G2_case3==2)*(M+1)))/P_G2_norm;
    %drop
    P_G2_drop(j)=(sum(P_G2(find(POLY_G2==0)))-P_G2_G3_drop(j)*P_G2_norm-P_G2_ack(j)*P_G2_norm)/P_G2_norm
    %check
    P_G2_check(j)=P_G2_ack(j)+P_G2_G2(j)+P_G2_G3_ret(j)+P_G2_G4(j)+ P_G2_G5_ret(j)+P_G2_G3_drop(j)+P_G2_G5_ts(j)+P_G2_G6_ret(j)+P_G2_drop(j)
    
    
    %%
%     %G3
%     P_G3=P_sta((M+1)*(M+2)+1:(M+1)*(M+3));
%     P_G3_norm=sum(P_G3);
%     P_G3_ack(j)=P_sta((M+1)*(M+3))/P_G3_norm;
%     POLY_G3=Policy_corre(j,(M+1)*(M+2)+1:(M+1)*(M+3));
%     P_G1_drop(j)=(sum(P_G1(find(POLY==0)))-P_G1(M+1))/P_G1_norm;
%     P_G1_ret(j)=sum(P_G1(find(POLY==2)))/P_G1_norm;
%     P_G1_ts(j)=sum(P_G1(find(POLY~=0&POLY~=2)))/P_G1_norm;
    
end

figure(1)

plot(SNR,P_G1_ack,'g',SNR,P_G1_drop,'b',SNR,P_G1_ret,'m',SNR,P_G1_ts,'k')

legend('P-G1-ack','P-G1-drop','P-G1-ret','P-G1-ts')

figure(2)

plot(SNR,P_G2_ack,'g',SNR,P_G2_G2,'b',SNR,P_G2_G3_ret,'m',SNR,P_G2_G4,'k',SNR,P_G2_G5_ret,'r',SNR,P_G2_G3_drop,'c',SNR,P_G2_G5_ts,'-.r'...
    ,SNR,P_G2_G6_ret,'-.c',SNR,P_G2_drop,'-.b')

legend('P-G2-ack','P-G2-G2','P-G2-G3-ret','P-G2-G4','P-G2-G5-ret','P-G2-G3-drop','P-G2-G5-ts','P-G2-G6-ret','P-G2-drop')
