clear;
clc;
%clf;


load AMC_policy_Perfect_Decoding_continuous_R_dop_0.075rho_0.8935.mat

SNRdB=[-5:30];
mu=0;
hI=1;
hR=1;
K=4;
N=1e7;
for ind_snr=1:length(SNRdB)
    
    ind_snr
    snr=10.^(SNRdB(ind_snr)./10);
    sigma=sqrt((1-rho).*snr./2);
    snr_m=snr_range(ind_snr,:);
    snr_m=[0 snr_m 1e100];
    policy_m=policy(ind_snr,:);
    policy_m=[policy_m(1) policy_m policy_m(end)];
    eR_m = normrnd(0,sigma,K,N);
    eI_m = normrnd(0,sigma,K,N);
    heR_m = normrnd(0,sqrt(snr./2),K,N);
    heI_m = normrnd(0,sqrt(snr./2),K,N);
    
    % K=1
    
    tt=1;
    eR = eR_m(tt,:);
    eI = eI_m(tt,:);
    heR = heR_m(tt,:);
    heI =heI_m(tt,:);
    gammae=heR.^2+ heI.^2;
    gamma=(sqrt(rho).*heR-eR).^2+(sqrt(rho).*heI-eI).^2;
    Optimal_rate1=interp1(snr_m,policy_m,gammae);
    Perfect_rate1=log2(1+gamma);
    needed_rate1=max(0,Optimal_rate1-Perfect_rate1);
    
    throughput_lharq(1,ind_snr)=sum(Optimal_rate1.*(Optimal_rate1<=Perfect_rate1))./N;
    throughput_harq(1,ind_snr)=throughput_lharq(1,ind_snr);
    outage(ind_snr)=sum((Optimal_rate1>Perfect_rate1))./N;
    
    %% K=2
    
    tt=2;
    eR = eR_m(tt,:);
    eI = eI_m(tt,:);
    heR = heR_m(tt,:);
    heI =heI_m(tt,:);
    gammae2=heR.^2+ heI.^2;
    gamma2=(sqrt(rho).*heR-eR).^2+(sqrt(rho).*heI-eI).^2;
    Optimal_rate2=interp1(snr_m,policy_m,gammae2);
    Perfect_rate2=log2(1+gamma2);
    needed_rate2=max(0,Optimal_rate2-Perfect_rate2);
    Second_T=(Optimal_rate1>Perfect_rate1).*(needed_rate1<=Optimal_rate2);
    throughput_lharq(2,ind_snr)=sum( Optimal_rate1.*(Optimal_rate1<=Perfect_rate1) + (Optimal_rate1>Perfect_rate1).*(needed_rate1<=Optimal_rate2).*...
        ( Optimal_rate2+(Optimal_rate1-needed_rate1).*(needed_rate1<=Optimal_rate2) ) .* ( Optimal_rate2<=Perfect_rate2) )...
        ./( N+sum( Second_T) );
    throughput_harq(2,ind_snr)=sum( Optimal_rate1.*(Optimal_rate1<=Perfect_rate1) + (Optimal_rate1>Perfect_rate1).*...
        Optimal_rate1.*(Optimal_rate1<=Perfect_rate1+Perfect_rate2) )...
        ./( N+sum( Optimal_rate1>Perfect_rate1 ) );
    
    %% K=3
    
    tt=3;
    eR = eR_m(tt,:);
    eI = eI_m(tt,:);
    heR = heR_m(tt,:);
    heI =heI_m(tt,:);
    gammae3=heR.^2+ heI.^2;
    gamma3=(sqrt(rho).*heR-eR).^2+(sqrt(rho).*heI-eI).^2;
    Optimal_rate3=interp1(snr_m,policy_m,gammae3);
    Perfect_rate3=log2(1+gamma3);
    needed_rate3=max(0,Optimal_rate3-Perfect_rate3);
    Third_T=(Optimal_rate1>Perfect_rate1).*(Optimal_rate2>Perfect_rate2).*(needed_rate1<=Optimal_rate2).*(needed_rate2<=Optimal_rate3);
    Third_T_C=(Optimal_rate1>Perfect_rate1).*(Optimal_rate1>Perfect_rate1+Perfect_rate2);
    
    throughput_lharq(3,ind_snr)=sum( Optimal_rate1.*(Optimal_rate1<=Perfect_rate1) + (Optimal_rate1>Perfect_rate1).*(needed_rate1<=Optimal_rate2).*...
        ( Optimal_rate2+(Optimal_rate1-needed_rate1).*(needed_rate1<=Optimal_rate2) ) .* ( Optimal_rate2<=Perfect_rate2) +...
        (Optimal_rate1>Perfect_rate1).*(Optimal_rate2>Perfect_rate2).*(needed_rate1<=Optimal_rate2).*(needed_rate2<=Optimal_rate3).*...
        ( Optimal_rate3+ ( Optimal_rate2-needed_rate2+(Optimal_rate1-needed_rate1).*(needed_rate1<=Optimal_rate2) ).*...
        (needed_rate2<=Optimal_rate3) ).*(Optimal_rate3<=Perfect_rate3) )  ...
        ./(N+sum(Second_T)+sum(Third_T));
    
    throughput_harq(3,ind_snr)=sum( Optimal_rate1.*(Optimal_rate1<=Perfect_rate1)+ (Optimal_rate1>Perfect_rate1).*...
        Optimal_rate1.*(Optimal_rate1<=Perfect_rate1+Perfect_rate2) +(Optimal_rate1>Perfect_rate1)...
        .*(Optimal_rate1>Perfect_rate1+Perfect_rate2).*Optimal_rate1.*(Optimal_rate1<=Perfect_rate1+Perfect_rate2+Perfect_rate3) )...
        ./( N+sum(Second_T)+sum(Third_T_C) );
    
     %% K=4
    
    tt=4;
    eR = eR_m(tt,:);
    eI = eI_m(tt,:);
    heR = heR_m(tt,:);
    heI =heI_m(tt,:);
    gammae4=heR.^2+ heI.^2;
    gamma4=(sqrt(rho).*heR-eR).^2+(sqrt(rho).*heI-eI).^2;
    Optimal_rate4=interp1(snr_m,policy_m,gammae4);
    Perfect_rate4=log2(1+gamma4);
    needed_rate4=max(0,Optimal_rate4-Perfect_rate4);
    Fourth_T=(Optimal_rate1>Perfect_rate1).*(Optimal_rate2>Perfect_rate2).*(Optimal_rate3>Perfect_rate3).*(needed_rate1<=Optimal_rate2)...
        .*(needed_rate2<=Optimal_rate3).*(needed_rate3<=Optimal_rate4);
    Fourth_T_C=(Optimal_rate1>Perfect_rate1).*(Optimal_rate1>Perfect_rate1+Perfect_rate2).*(Optimal_rate1>Perfect_rate1+Perfect_rate2+Perfect_rate3);
    
    throughput_lharq(4,ind_snr)=sum( Optimal_rate1.*(Optimal_rate1<=Perfect_rate1) + (Optimal_rate1>Perfect_rate1).*(needed_rate1<=Optimal_rate2).*...
        ( Optimal_rate2+(Optimal_rate1-needed_rate1).*(needed_rate1<=Optimal_rate2) ) .* ( Optimal_rate2<=Perfect_rate2) +...
        (Optimal_rate1>Perfect_rate1).*(Optimal_rate2>Perfect_rate2).*(needed_rate1<=Optimal_rate2).*(needed_rate2<=Optimal_rate3).*...
        ( Optimal_rate3+ ( Optimal_rate2-needed_rate2+(Optimal_rate1-needed_rate1).*(needed_rate1<=Optimal_rate2) ).*...
        (needed_rate2<=Optimal_rate3) ).*(Optimal_rate3<=Perfect_rate3)+  ...
        (Optimal_rate1>Perfect_rate1).*(Optimal_rate2>Perfect_rate2).*(Optimal_rate3>Perfect_rate3).*(needed_rate1<=Optimal_rate2).*(needed_rate2<=Optimal_rate3).*(needed_rate3<=Optimal_rate4).*...
        ( Optimal_rate4+ ( Optimal_rate3-needed_rate3+( Optimal_rate2-needed_rate2+(Optimal_rate1-needed_rate1).*(needed_rate1<=Optimal_rate2) ).*...
        (needed_rate2<=Optimal_rate3) ).*(needed_rate3<=Optimal_rate4) ).*(Optimal_rate4<=Perfect_rate4))...
        ./(N+sum(Second_T)+sum(Third_T)+sum(Fourth_T));
    
    throughput_harq(4,ind_snr)=sum( Optimal_rate1.*(Optimal_rate1<=Perfect_rate1)+ (Optimal_rate1>Perfect_rate1).*...
        Optimal_rate1.*(Optimal_rate1<=Perfect_rate1+Perfect_rate2) +(Optimal_rate1>Perfect_rate1)...
        .*(Optimal_rate1>Perfect_rate1+Perfect_rate2).*Optimal_rate1.*(Optimal_rate1<=Perfect_rate1+Perfect_rate2+Perfect_rate3)...
        +(Optimal_rate1>Perfect_rate1).*(Optimal_rate1>Perfect_rate1+Perfect_rate2).* (Optimal_rate1>Perfect_rate1+Perfect_rate2+Perfect_rate3).*...
        Optimal_rate1.*(Optimal_rate1<=Perfect_rate1+Perfect_rate2+Perfect_rate3+Perfect_rate4) )...
        ./( N+sum(Second_T)+sum(Third_T_C)+sum(Fourth_T_C) );
    
    
end

save(['LHARQ_AMC_policy_Perfect_Decoding_continuous_R_dop_',num2str(dop),'_rho_',num2str(rho),'.mat'])

% figure(1)
% plot(SNRdB,throughput_lharq(1,:),'r',SNRdB,throughput_lharq(2,:),'g',SNRdB,throughput_lharq(3,:),'m',SNRdB,throughput_lharq(4,:),'b')
% hold on
% plot(SNRdB,throughput_harq(2,:),'-.b',SNRdB,throughput_harq(3,:),'-.m',SNRdB,throughput_harq(4,:),'-.g')

% hold on
% plot(SNRdB,throughput_HARQ(2,:),'-.r',SNRdB,throughput_HARQ(3,:),'-.b',SNRdB,ergodic_capacity,'-.m')
% 
% legend('AMC','LHARQ, K=2','LHARQ, K=3','HARQ, K=2','HARQ, K=3', 'Ergodic capacity')
% xlabel('SNR [dB]')
% ylabel('Throughput')
% title(['rho=',num2str(rho),])
% grid
% 
% figure(2)
% semilogy(SNRdB,outage)
% 
% xlabel('SNR [dB]')
% ylabel('Outage')
% title(['rho=',num2str(rho),])
% grid


% snr_m=snr_range(ind_snr,:);
%     snr_m=[0 snr_m 1e100];
%     policy_m=policy(ind_snr,:);
%     policy_m=[policy_m(1) policy_m policy_m(end)];
%     eR_m = normrnd(0,sigma,K,N);
%     eI_m = normrnd(0,sigma,K,N);
%     heR_m = normrnd(0,sqrt(snr./2),K,N);
%     heI_m = normrnd(0,sqrt(snr./2),K,N);
%     gammae=heR_m.^2+ heI_m.^2;
%     gamma=(sqrt(rho).*heR_m-eR_m).^2+(sqrt(rho).*heI_m-eI_m).^2;
%     gammath=log2(1+gamma);
% 
%     Optimal_rate=interp1(snr_m,policy_m,gammae);
%     bb=Optimal_rate(1,:);
%     Perfect_rate= gammath;
%     temp=Optimal_rate-Perfect_rate;
%     Needed_rate=zeros(K,N);
%     Needed_rate(2:K,1:N)=temp(1:K-1,1:N);
%     Needed_rate=max(0,Needed_rate);
%     %% L-HARQ
%     Plus_round=Optimal_rate>gammath;
%     Success_round=~Plus_round;
%     Plus_round= cumprod(Plus_round);
%     exist_T=ones(K,N);
%     exist_T(2:K,1:N)=Plus_round(1:K-1,1:N);
%     Number_exp=sum(Plus_round,[2]);
%     Number_exp=[N;Number_exp(1:K-1)];
%     Number_exp=cumsum(Number_exp);
%     % reward after direct decoding
%    
%     Reward_d=Optimal_rate.*Success_round;
%     Reward_d=Reward_d.*exist_T;
%     Reward_d=cumsum(Reward_d);
%     % reward of backtrack decoding
%     gammath_b=zeros(K,N);
%     gammath_b(2:K,1:N)=gammath(1:K-1,1:N);
%     (Needed_rate<=Optimal_rate)
%     cumprod((Needed_rate<=Optimal_rate),'reverse')
%     exist_T
%     
%     Reward_b=gammath_b.*exist_T.*(Needed_rate<=Optimal_rate);
%     Success_past=cumsum(Success_round);
%     Reward_b=cumsum(Reward_b).*(Success_past>0);
%     throughput_lharq(:,ind_snr)=sum(Reward_d+Reward_b,[2])./Number_exp;
%     %% HARQ
