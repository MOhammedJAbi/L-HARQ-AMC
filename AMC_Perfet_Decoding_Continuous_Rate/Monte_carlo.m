clear;
clc;
clf;

%load AMC_Perfect_Decoding_continuous_R_rho_0.99.mat
load AMC_Perfect_Decoding_continuous_R_rho_0.9V2.mat

SNRdB=[-5:30];
mu=0;
hI=1;
hR=1;
K=3;
N=1e6;
for ind_snr=1:length(SNRdB)
    ind_snr
    snr=10.^(SNRdB(ind_snr)./10);
    sigma=sqrt((1-rho).*snr./2);
    %% K=1
    eR = normrnd(0,sigma,1,N);
    eI = normrnd(0,sigma,1,N);
    heR = normrnd(0,sqrt(snr./2),1,N);
    heI = normrnd(0,sqrt(snr./2),1,N);
    gammae=heR.^2+ heI.^2;
    gamma=(sqrt(rho).*heR-eR).^2+(sqrt(rho).*heI-eI).^2;
    snr_m=snr_range(ind_snr,:);
    snr_m=[0 snr_m 1e100];
    policy_m=policy(ind_snr,:);
    policy_m=[policy_m(1) policy_m policy_m(end)];
    Optimal_rate=interp1(snr_m,policy_m,gammae);
    Perfect_rate=interp1(snr_m,policy_m,gamma);
    needed_rate=max(0,Optimal_rate-Perfect_rate);

    throughput_LHARQ(1,ind_snr)=sum(Optimal_rate.*(Optimal_rate<=log2(gamma+1)))./N;
    throughput_HARQ(1,ind_snr)=throughput_LHARQ(1,ind_snr);
    outage(ind_snr)=sum((Optimal_rate>log2(gamma+1)))./N;
    
     %% K=2
    
    eR = normrnd(0,sigma,1,N);
    eI = normrnd(0,sigma,1,N);
    heR = normrnd(0,sqrt(snr./2),1,N);
    heI = normrnd(0,sqrt(snr./2),1,N);
    gammae2=heR.^2+ heI.^2;
    gamma2=(sqrt(rho).*heR-eR).^2+(sqrt(rho).*heI-eI).^2;
    Optimal_rate2=interp1(snr_m,policy_m,gammae2);
    Perfect_rate2=interp1(snr_m,policy_m,gamma);
    needed_rate2=max(0,Optimal_rate2-Perfect_rate2);
    
    throughput_LHARQ(2,ind_snr)=sum(Optimal_rate.*(Optimal_rate<=log2(gamma+1))+(Optimal_rate>log2(gamma+1)).*...
        (Optimal_rate2+(Optimal_rate-needed_rate).*(needed_rate<=Optimal_rate2)).*(Optimal_rate2<=log2(gamma2+1)) )...
        ./(N+sum(Optimal_rate>log2(gamma+1)));
    throughput_HARQ(2,ind_snr)=sum(Optimal_rate.*(Optimal_rate<=log2(gamma+1))+(Optimal_rate>log2(gamma+1)).*...
        Optimal_rate.*(Optimal_rate<=log2(gamma2+1)+log2(gamma+1)) )...
        ./(N+sum(Optimal_rate>log2(gamma+1)));
    
    %% K=3
    
    eR = normrnd(0,sigma,1,N);
    eI = normrnd(0,sigma,1,N);
    heR = normrnd(0,sqrt(snr./2),1,N);
    heI = normrnd(0,sqrt(snr./2),1,N);
    gammae3=heR.^2+ heI.^2;
    gamma3=(sqrt(rho).*heR-eR).^2+(sqrt(rho).*heI-eI).^2;
    Optimal_rate3=interp1(snr_m,policy_m,gammae3);
    Second_T=Optimal_rate>log2(gamma+1);
    Third_T=(Optimal_rate>log2(gamma+1)).*(Optimal_rate2>log2(gamma2+1));
    Third_T_C=(Optimal_rate>log2(gamma+1)).*(Optimal_rate>log2(gamma2+1)+log2(gamma+1));
    throughput_LHARQ(3,ind_snr)=sum(Optimal_rate.*(Optimal_rate<=log2(gamma+1))+(Optimal_rate>log2(gamma+1)).*...
        (Optimal_rate2+(Optimal_rate-needed_rate).*(needed_rate<=Optimal_rate2)).*(Optimal_rate2<=log2(gamma2+1))+...
        (Optimal_rate>log2(gamma+1)).*(Optimal_rate2>log2(gamma2+1)).*(Optimal_rate3+(Optimal_rate2-needed_rate2+...
        (Optimal_rate-needed_rate).*(needed_rate<=Optimal_rate2)).*(needed_rate2<=Optimal_rate3)).*(Optimal_rate3<=log2(gamma3+1)))...
        ./(N+sum(Second_T)+sum(Third_T));
    throughput_HARQ(3,ind_snr)=sum(Optimal_rate.*(Optimal_rate<=log2(gamma+1))+(Optimal_rate>log2(gamma+1)).*...
        Optimal_rate.*(Optimal_rate<=log2(gamma2+1)+log2(gamma+1)) +(Optimal_rate>log2(gamma+1))...
        .*(Optimal_rate>log2(gamma2+1)+log2(gamma+1)).*Optimal_rate.*(Optimal_rate<=log2(gamma2+1)+log2(gamma+1)+log2(gamma3+1)))...
        ./(N+sum(Optimal_rate>log2(gamma+1))+sum(Third_T_C));
    
    
end

save(['LHARQ_Perfect_Decoding_continuous_R_rho_',num2str(rho),'.mat'])

figure(1)
plot(SNRdB,throughput_LHARQ(1,:),'r',SNRdB,throughput_LHARQ(2,:),'b',SNRdB,throughput_LHARQ(3,:),'g')
hold on
plot(SNRdB,throughput_HARQ(2,:),'-.r',SNRdB,throughput_HARQ(3,:),'-.b',SNRdB,ergodic_capacity,'-.m')

legend('AMC','LHARQ, K=2','LHARQ, K=3','HARQ, K=2','HARQ, K=3', 'Ergodic capacity')
xlabel('SNR [dB]')
ylabel('Throughput')
title(['rho=',num2str(rho),])
grid
% 
% figure(2)
% semilogy(SNRdB,outage)
% 
% xlabel('SNR [dB]')
% ylabel('Outage')
% title(['rho=',num2str(rho),])
% grid
