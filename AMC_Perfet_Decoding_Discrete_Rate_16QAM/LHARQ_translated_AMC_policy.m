clear;
clc;
%clf;


load AMC_policy_16QAM_Perfect_Decoding_continuous_R_dop_0.05rho_0.95156.mat

load MI_16QAM.mat
MI_16QAM=[0 MI_16QAM(1:end-1)' 4];
th=[0 10.^((SNR(1:end-1))'./10) 1e10];
Qam_th=@(R) interp1(MI_16QAM,th,R);
Qam_r=@(gamma) interp1(th,MI_16QAM,gamma);

SNRdB=[-5:30];
mu=0;
hI=1;
hR=1;
K=4;
N=1e6;

for ind_snr=1:length(SNRdB)
    
    ind_snr
    snr=10.^(SNRdB(ind_snr)./10);
    sigma=sqrt((1-rho).*snr./2);
    snr_m=snr_range(ind_snr,:);
    %snr_m=[0 snr_m 1e100];
    policy_m=policy(ind_snr,:);
    %policy_m=[policy_m(1) policy_m policy_m(end)];
    policy_m=[policy_m policy_m(end)];
    eR_m = normrnd(0,sigma,K,N);
    eI_m = normrnd(0,sigma,K,N);
    heR_m = normrnd(0,sqrt(snr./2),K,N);
    heI_m = normrnd(0,sqrt(snr./2),K,N);
    delta_set=[0 logspace(log10(snr_range(ind_snr,1)),log10(snr_range(ind_snr,end)),50)];
    
   
    
    for ind_delta=1:length(delta_set)

        snr_mm=snr_m-delta_set(ind_delta);
       snr_mm=[snr_mm 1e100];
         % K=1
        tt=1;
        eR = eR_m(tt,:);
        eI = eI_m(tt,:);
        heR = heR_m(tt,:);
        heI =heI_m(tt,:);
        gammae=heR.^2+ heI.^2;
        gamma=(sqrt(rho).*heR-eR).^2+(sqrt(rho).*heI-eI).^2;
        Optimal_rate1=interp1(snr_mm,policy_m,gammae);
        Perfect_rate1=Qam_r(gamma);%log2(1+gamma);
        needed_rate1=max(0,Optimal_rate1-Perfect_rate1);
        
        throughput_lharq_temp(1,ind_delta)=sum(Optimal_rate1.*(Optimal_rate1<=Perfect_rate1))./N;
        throughput_harq_temp(1,ind_delta)=throughput_lharq_temp(1,ind_delta);
        outage_temp(ind_delta)=sum((Optimal_rate1>Perfect_rate1))./N;
        
        %% K=2
        
        tt=2;
        eR = eR_m(tt,:);
        eI = eI_m(tt,:);
        heR = heR_m(tt,:);
        heI =heI_m(tt,:);
        gammae2=heR.^2+ heI.^2;
        gamma2=(sqrt(rho).*heR-eR).^2+(sqrt(rho).*heI-eI).^2;
        Optimal_rate2=interp1(snr_mm,policy_m,gammae2);
        Perfect_rate2=Qam_r(gamma2);%log2(1+gamma2);
        needed_rate2=max(0,Optimal_rate2-Perfect_rate2);
        Second_T=(Optimal_rate1>Perfect_rate1).*(needed_rate1<=Optimal_rate2);
        throughput_lharq_temp(2,ind_delta)=sum( Optimal_rate1.*(Optimal_rate1<=Perfect_rate1) + (Optimal_rate1>Perfect_rate1).*(needed_rate1<=Optimal_rate2).*...
            ( Optimal_rate2+(Optimal_rate1-needed_rate1).*(needed_rate1<=Optimal_rate2) ) .* ( Optimal_rate2<=Perfect_rate2) )...
            ./( N+sum( Second_T) );
        throughput_harq_temp(2,ind_delta)=sum( Optimal_rate1.*(Optimal_rate1<=Perfect_rate1) + (Optimal_rate1>Perfect_rate1).*...
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
        Optimal_rate3=interp1(snr_mm,policy_m,gammae3);
        Perfect_rate3=Qam_r(gamma3); %log2(1+gamma3);
        needed_rate3=max(0,Optimal_rate3-Perfect_rate3);
        Third_T=(Optimal_rate1>Perfect_rate1).*(Optimal_rate2>Perfect_rate2).*(needed_rate1<=Optimal_rate2).*(needed_rate2<=Optimal_rate3);
        Third_T_C=(Optimal_rate1>Perfect_rate1).*(Optimal_rate1>Perfect_rate1+Perfect_rate2);
        
        throughput_lharq_temp(3,ind_delta)=sum( Optimal_rate1.*(Optimal_rate1<=Perfect_rate1) + (Optimal_rate1>Perfect_rate1).*(needed_rate1<=Optimal_rate2).*...
            ( Optimal_rate2+(Optimal_rate1-needed_rate1).*(needed_rate1<=Optimal_rate2) ) .* ( Optimal_rate2<=Perfect_rate2) +...
            (Optimal_rate1>Perfect_rate1).*(Optimal_rate2>Perfect_rate2).*(needed_rate1<=Optimal_rate2).*(needed_rate2<=Optimal_rate3).*...
            ( Optimal_rate3+ ( Optimal_rate2-needed_rate2+(Optimal_rate1-needed_rate1).*(needed_rate1<=Optimal_rate2) ).*...
            (needed_rate2<=Optimal_rate3) ).*(Optimal_rate3<=Perfect_rate3) )  ...
            ./(N+sum(Second_T)+sum(Third_T));
        
        throughput_harq_temp(3,ind_delta)=sum( Optimal_rate1.*(Optimal_rate1<=Perfect_rate1)+ (Optimal_rate1>Perfect_rate1).*...
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
        Optimal_rate4=interp1(snr_mm,policy_m,gammae4);
        Perfect_rate4=Qam_r(gamma4); %log2(1+gamma4);
        needed_rate4=max(0,Optimal_rate4-Perfect_rate4);
        Fourth_T=(Optimal_rate1>Perfect_rate1).*(Optimal_rate2>Perfect_rate2).*(Optimal_rate3>Perfect_rate3).*(needed_rate1<=Optimal_rate2)...
            .*(needed_rate2<=Optimal_rate3).*(needed_rate3<=Optimal_rate4);
        Fourth_T_C=(Optimal_rate1>Perfect_rate1).*(Optimal_rate1>Perfect_rate1+Perfect_rate2).*(Optimal_rate1>Perfect_rate1+Perfect_rate2+Perfect_rate3);
        
        throughput_lharq_temp(4,ind_delta)=sum( Optimal_rate1.*(Optimal_rate1<=Perfect_rate1) + (Optimal_rate1>Perfect_rate1).*(needed_rate1<=Optimal_rate2).*...
            ( Optimal_rate2+(Optimal_rate1-needed_rate1).*(needed_rate1<=Optimal_rate2) ) .* ( Optimal_rate2<=Perfect_rate2) +...
            (Optimal_rate1>Perfect_rate1).*(Optimal_rate2>Perfect_rate2).*(needed_rate1<=Optimal_rate2).*(needed_rate2<=Optimal_rate3).*...
            ( Optimal_rate3+ ( Optimal_rate2-needed_rate2+(Optimal_rate1-needed_rate1).*(needed_rate1<=Optimal_rate2) ).*...
            (needed_rate2<=Optimal_rate3) ).*(Optimal_rate3<=Perfect_rate3)+  ...
            (Optimal_rate1>Perfect_rate1).*(Optimal_rate2>Perfect_rate2).*(Optimal_rate3>Perfect_rate3).*(needed_rate1<=Optimal_rate2).*(needed_rate2<=Optimal_rate3).*(needed_rate3<=Optimal_rate4).*...
            ( Optimal_rate4+ ( Optimal_rate3-needed_rate3+( Optimal_rate2-needed_rate2+(Optimal_rate1-needed_rate1).*(needed_rate1<=Optimal_rate2) ).*...
            (needed_rate2<=Optimal_rate3) ).*(needed_rate3<=Optimal_rate4) ).*(Optimal_rate4<=Perfect_rate4))...
            ./(N+sum(Second_T)+sum(Third_T)+sum(Fourth_T));
        
        throughput_harq_temp(4,ind_delta)=sum( Optimal_rate1.*(Optimal_rate1<=Perfect_rate1)+ (Optimal_rate1>Perfect_rate1).*...
            Optimal_rate1.*(Optimal_rate1<=Perfect_rate1+Perfect_rate2) +(Optimal_rate1>Perfect_rate1)...
            .*(Optimal_rate1>Perfect_rate1+Perfect_rate2).*Optimal_rate1.*(Optimal_rate1<=Perfect_rate1+Perfect_rate2+Perfect_rate3)...
            +(Optimal_rate1>Perfect_rate1).*(Optimal_rate1>Perfect_rate1+Perfect_rate2).* (Optimal_rate1>Perfect_rate1+Perfect_rate2+Perfect_rate3).*...
            Optimal_rate1.*(Optimal_rate1<=Perfect_rate1+Perfect_rate2+Perfect_rate3+Perfect_rate4) )...
            ./( N+sum(Second_T)+sum(Third_T_C)+sum(Fourth_T_C) );
    end
    
    for ind_k=1:K
        [throughput_lharq(ind_k,ind_snr),ind] = max(throughput_lharq_temp(ind_k,:));
        delta_lharq(ind_k,ind_snr)=delta_set(ind);
        outage(ind_k,ind_snr)=outage_temp(ind);
        [throughput_harq(ind_k,ind_snr),ind] = max(throughput_harq_temp(ind_k,:));
        delta_harq(ind_k,ind_snr)=delta_set(ind);
    end
    
    
end
save LHARQ_translated_AMC_policy_16QAM_Perfect_Decoding_continuous_R_dop_0.075rho_0.9_V_leger.mat  dop rho SNRdB ...
    throughput_lharq throughput_harq outage snr_range policy delta_harq delta_lharq

