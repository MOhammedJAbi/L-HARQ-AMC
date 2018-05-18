clear;
clc;
%clf;


load AMC_policy_16_64QAM_Perfect_Decoding_continuous_R_dop_0.1rho_0.8167.mat

load MI_64QAM_SNRdB=-20_30_sigmap2=0_0.mat
MI_64QAM=[0 MI_PN_BLT(1:end-5)' 6];
th_64=[0 10.^((SNRdB(1:end-5))./10) 1e10];
Qam_th_64=@(R) interp1(MI_64QAM,th_64,R)
Qam_r_64=@(gamma) interp1(th_64,MI_64QAM,gamma)

load MI_16QAM.mat
MI_16QAM=[0 MI_16QAM(1:end-1)' 4];
th_16=[0 10.^((SNR(1:end-1))'./10) 1e10];
Qam_th_16=@(R) interp1(MI_16QAM,th_16,R)
Qam_r_16=@(gamma) interp1(th_16,MI_16QAM,gamma)


SNRdB=[-5:40];
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
    snr_m=[0 snr_m 1e100];
    policy_m=policy(ind_snr,:);
    modulation_m=modulation(ind_snr,:);
    modulation_m=[modulation_m(1) modulation_m modulation_m(end)];
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
    Optimal_modulation1=interp1(snr_m,modulation_m,gammae,'nearest');
    for ii=1:N
        if Optimal_modulation1(ii)==1
            Perfect_rate1(ii)=Qam_r_16(gamma(ii));%log2(1+gamma);
        elseif Optimal_modulation1(ii)==2
            Perfect_rate1(ii)=Qam_r_64(gamma(ii));
        else
            error('modulation not supported')
        end
    end
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
    Optimal_modulation2=interp1(snr_m,modulation_m,gammae2,'nearest');
    for ii=1:N
        if Optimal_modulation2(ii)==1
            Perfect_rate2(ii)=Qam_r_16(gamma2(ii));%log2(1+gamma);
        elseif Optimal_modulation2(ii)==2
            Perfect_rate2(ii)=Qam_r_64(gamma2(ii));
        else
            error('modulation not supported')
        end
    end
    %Perfect_rate2=Qam_r(gamma2);%log2(1+gamma2);
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
    Optimal_modulation3=interp1(snr_m,modulation_m,gammae3,'nearest');
    for ii=1:N
        if Optimal_modulation3(ii)==1
            Perfect_rate3(ii)=Qam_r_16(gamma3(ii));%log2(1+gamma);
        elseif Optimal_modulation3(ii)==2
            Perfect_rate3(ii)=Qam_r_64(gamma3(ii));
        else
            error('modulation not supported')
        end
    end
    %Perfect_rate3=Qam_r(gamma3); %log2(1+gamma3);
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
    Optimal_modulation4=interp1(snr_m,modulation_m,gammae4,'nearest');
    for ii=1:N
        if Optimal_modulation4(ii)==1
            Perfect_rate4(ii)=Qam_r_16(gamma4(ii));%log2(1+gamma);
        elseif Optimal_modulation4(ii)==2
            Perfect_rate4(ii)=Qam_r_64(gamma4(ii));
        else
            error('modulation not supported')
        end
    end
    %Perfect_rate4=Qam_r(gamma4); %log2(1+gamma4);
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
save LHARQ_AMC_policy_16_64QAM_Perfect_Decoding_continuous_R_dop_0.1rho_0.8_V_leger.mat  dop rho SNRdB ...
    throughput_lharq throughput_harq outage snr_range policy modulation

