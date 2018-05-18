clear;
%clc;
%clf;

load AMC_policy_64QAM_Perfect_Decoding_Discrete_R_dop_0.1rho_0.8167.mat

load MI_64QAM_SNRdB=-20_30_sigmap2=0_0.mat
MI_16QAM=[0 MI_PN_BLT(1:end-5)' 6];
th=[0 10.^((SNRdB(1:end-5))./10) 1e10];
Qam_th=@(R) interp1(MI_16QAM,th,R);
Qam_r=@(gamma) interp1(th,MI_16QAM,gamma);

SNRdB=[-5:40];
K=2;
N_realization=1e5;

for ind_snr=1:1:length(SNRdB)
    
    snr=10.^(SNRdB(ind_snr)./10);
    sigma=sqrt((1-rho).*snr./2);
    snr_m=[1e-20 snr_range(ind_snr,:) 1e100];
    policy_m=[policy(ind_snr,1) policy(ind_snr,:) policy(ind_snr,end)];
    eR = normrnd(0,sigma,1,N_realization);
    eI = normrnd(0,sigma,1,N_realization);
    heR = normrnd(0,sqrt(snr./2),1,N_realization);
    heI = normrnd(0,sqrt(snr./2),1,N_realization);
    gammae=heR.^2+ heI.^2;
    gamma=(sqrt(rho).*heR-eR).^2+(sqrt(rho).*heI-eI).^2;
    
    delta_set=logspace(0,log10(snr_range(ind_snr,end-1)),50);
    
    for ind_delta=1:length(delta_set)
        [ind_snr ind_delta]
        snr_mm=snr_m./delta_set(ind_delta);
        counter=0;
        reward=0;
        while counter<N_realization-1
            counter=counter+1;
            rate_round_1=interp1(snr_mm,policy_m,gammae(counter),'nearest');
            wep_round_1=Qam_r(gamma(counter))<rate_round_1;
            ack_round_1=1-wep_round_1;
            if ack_round_1==1
                reward=reward+rate_round_1;
            else
                BK_rate_round_1 = rate_round_1-Qam_r(gamma(counter));
                rate_round_2=interp1(snr_mm,policy_m,gammae(counter+1),'nearest');
                if BK_rate_round_1<=rate_round_2
                    wep_round_2=Qam_r(gamma(counter+1))<rate_round_2;
                    ack_round_2=1-wep_round_2;
                    reward=reward+(rate_round_2+rate_round_1-BK_rate_round_1)*ack_round_2;
                    counter=counter+1;
                end
            end
        end
        throughput_lharq_delta(ind_delta,ind_snr)=reward/counter;
    end
    [throughput_lharq(ind_snr),jj]=max(throughput_lharq_delta(:,ind_snr));
    optimal_delta(ind_snr)=delta_set(jj);
end


save LHARQ_translated_AMC_policy_64QAM_K_2_Perfect_Decoding_Discrete_R_dop_0.1_rho_0.8_V_LHARQ_Correct.mat  dop rho SNRdB ...
    throughput_lharq_delta throughput_lharq  snr_range policy optimal_delta

