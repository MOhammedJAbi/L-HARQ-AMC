clear;
%clc;
%clf;


[FInterpolation1,FInterpolation2,FInterpolation3,FInterpolation4] = WEP_space(1);

load AMC_policy_16QAM_Turbo_code_dop_0.05rho_0.95156.mat

SNRdB=[-5:30];
K=1;
N_realization=1.25e5;

for ind_snr=1:length(SNRdB)
    ind_snr
    snr=10.^(SNRdB(ind_snr)./10);
    sigma=sqrt((1-rho).*snr./2);
    snr_m=snr_range(ind_snr,:);
    policy_m=policy(ind_snr,:);
    eR = normrnd(0,sigma,1,N_realization);
    eI = normrnd(0,sigma,1,N_realization);
    heR = normrnd(0,sqrt(snr./2),1,N_realization);
    heI = normrnd(0,sqrt(snr./2),1,N_realization);
    gammae=heR.^2+ heI.^2;
    gamma=(sqrt(rho).*heR-eR).^2+(sqrt(rho).*heI-eI).^2;
    
    counter=0;
    reward=0;
    number_ack_1=0;
    while counter<N_realization
        counter=counter+1;
        rate_round_1=interp1(snr_m,policy_m,gammae(counter),'nearest');
        wep_round_1=WEP_function_Fast(rate_round_1,0,gamma(counter),FInterpolation1,FInterpolation2,FInterpolation3,FInterpolation4);
        ack_round_1=binornd(1, 1-wep_round_1);
        if ack_round_1==1
            reward=reward+rate_round_1;
        end
    end
    throughput_amc(ind_snr)= reward/counter;
end


save LHARQ_translated_AMC_policy_16QAM_turbo_code_K_1_dop_0.05_rho_0.95.mat  dop rho SNRdB throughput_amc



%         Perfect_rate1=Qam_r(gamma);%log2(1+gamma);
%         needed_rate1=max(0,Optimal_rate1-Perfect_rate1);
%
%         throughput_lharq_temp(1,ind_delta)=sum(Optimal_rate1.*(Optimal_rate1<=Perfect_rate1))./N;
%         throughput_harq_temp(1,ind_delta)=throughput_lharq_temp(1,ind_delta);
%         outage_temp(ind_delta)=sum((Optimal_rate1>Perfect_rate1))./N;
%
%         %% K=2
%
%         tt=2;
%         eR = eR_m(tt,:);
%         eI = eI_m(tt,:);
%         heR = heR_m(tt,:);
%         heI =heI_m(tt,:);
%         gammae2=heR.^2+ heI.^2;
%         gamma2=(sqrt(rho).*heR-eR).^2+(sqrt(rho).*heI-eI).^2;
%         Optimal_rate2=interp1(snr_mm,policy_m,gammae2);
%         Perfect_rate2=Qam_r(gamma2);%log2(1+gamma2);
%         needed_rate2=max(0,Optimal_rate2-Perfect_rate2);
%         Second_T=(Optimal_rate1>Perfect_rate1).*(needed_rate1<=Optimal_rate2);
%         throughput_lharq_temp(2,ind_delta)=sum( Optimal_rate1.*(Optimal_rate1<=Perfect_rate1) + (Optimal_rate1>Perfect_rate1).*(needed_rate1<=Optimal_rate2).*...
%             ( Optimal_rate2+(Optimal_rate1-needed_rate1).*(needed_rate1<=Optimal_rate2) ) .* ( Optimal_rate2<=Perfect_rate2) )...
%             ./( N+sum( Second_T) );
%         throughput_harq_temp(2,ind_delta)=sum( Optimal_rate1.*(Optimal_rate1<=Perfect_rate1) + (Optimal_rate1>Perfect_rate1).*...
%             Optimal_rate1.*(Optimal_rate1<=Perfect_rate1+Perfect_rate2) )...
%             ./( N+sum( Optimal_rate1>Perfect_rate1 ) );
%
%         %% K=3
%
%         tt=3;
%         eR = eR_m(tt,:);
%         eI = eI_m(tt,:);
%         heR = heR_m(tt,:);
%         heI =heI_m(tt,:);
%         gammae3=heR.^2+ heI.^2;
%         gamma3=(sqrt(rho).*heR-eR).^2+(sqrt(rho).*heI-eI).^2;
%         Optimal_rate3=interp1(snr_mm,policy_m,gammae3);
%         Perfect_rate3=Qam_r(gamma3); %log2(1+gamma3);
%         needed_rate3=max(0,Optimal_rate3-Perfect_rate3);
%         Third_T=(Optimal_rate1>Perfect_rate1).*(Optimal_rate2>Perfect_rate2).*(needed_rate1<=Optimal_rate2).*(needed_rate2<=Optimal_rate3);
%         Third_T_C=(Optimal_rate1>Perfect_rate1).*(Optimal_rate1>Perfect_rate1+Perfect_rate2);
%
%         throughput_lharq_temp(3,ind_delta)=sum( Optimal_rate1.*(Optimal_rate1<=Perfect_rate1) + (Optimal_rate1>Perfect_rate1).*(needed_rate1<=Optimal_rate2).*...
%             ( Optimal_rate2+(Optimal_rate1-needed_rate1).*(needed_rate1<=Optimal_rate2) ) .* ( Optimal_rate2<=Perfect_rate2) +...
%             (Optimal_rate1>Perfect_rate1).*(Optimal_rate2>Perfect_rate2).*(needed_rate1<=Optimal_rate2).*(needed_rate2<=Optimal_rate3).*...
%             ( Optimal_rate3+ ( Optimal_rate2-needed_rate2+(Optimal_rate1-needed_rate1).*(needed_rate1<=Optimal_rate2) ).*...
%             (needed_rate2<=Optimal_rate3) ).*(Optimal_rate3<=Perfect_rate3) )  ...
%             ./(N+sum(Second_T)+sum(Third_T));
%
%         throughput_harq_temp(3,ind_delta)=sum( Optimal_rate1.*(Optimal_rate1<=Perfect_rate1)+ (Optimal_rate1>Perfect_rate1).*...
%             Optimal_rate1.*(Optimal_rate1<=Perfect_rate1+Perfect_rate2) +(Optimal_rate1>Perfect_rate1)...
%             .*(Optimal_rate1>Perfect_rate1+Perfect_rate2).*Optimal_rate1.*(Optimal_rate1<=Perfect_rate1+Perfect_rate2+Perfect_rate3) )...
%             ./( N+sum(Second_T)+sum(Third_T_C) );
%
%         %% K=4
%
%         tt=4;
%         eR = eR_m(tt,:);
%         eI = eI_m(tt,:);
%         heR = heR_m(tt,:);
%         heI =heI_m(tt,:);
%         gammae4=heR.^2+ heI.^2;
%         gamma4=(sqrt(rho).*heR-eR).^2+(sqrt(rho).*heI-eI).^2;
%         Optimal_rate4=interp1(snr_mm,policy_m,gammae4);
%         Perfect_rate4=Qam_r(gamma4); %log2(1+gamma4);
%         needed_rate4=max(0,Optimal_rate4-Perfect_rate4);
%         Fourth_T=(Optimal_rate1>Perfect_rate1).*(Optimal_rate2>Perfect_rate2).*(Optimal_rate3>Perfect_rate3).*(needed_rate1<=Optimal_rate2)...
%             .*(needed_rate2<=Optimal_rate3).*(needed_rate3<=Optimal_rate4);
%         Fourth_T_C=(Optimal_rate1>Perfect_rate1).*(Optimal_rate1>Perfect_rate1+Perfect_rate2).*(Optimal_rate1>Perfect_rate1+Perfect_rate2+Perfect_rate3);
%
%         throughput_lharq_temp(4,ind_delta)=sum( Optimal_rate1.*(Optimal_rate1<=Perfect_rate1) + (Optimal_rate1>Perfect_rate1).*(needed_rate1<=Optimal_rate2).*...
%             ( Optimal_rate2+(Optimal_rate1-needed_rate1).*(needed_rate1<=Optimal_rate2) ) .* ( Optimal_rate2<=Perfect_rate2) +...
%             (Optimal_rate1>Perfect_rate1).*(Optimal_rate2>Perfect_rate2).*(needed_rate1<=Optimal_rate2).*(needed_rate2<=Optimal_rate3).*...
%             ( Optimal_rate3+ ( Optimal_rate2-needed_rate2+(Optimal_rate1-needed_rate1).*(needed_rate1<=Optimal_rate2) ).*...
%             (needed_rate2<=Optimal_rate3) ).*(Optimal_rate3<=Perfect_rate3)+  ...
%             (Optimal_rate1>Perfect_rate1).*(Optimal_rate2>Perfect_rate2).*(Optimal_rate3>Perfect_rate3).*(needed_rate1<=Optimal_rate2).*(needed_rate2<=Optimal_rate3).*(needed_rate3<=Optimal_rate4).*...
%             ( Optimal_rate4+ ( Optimal_rate3-needed_rate3+( Optimal_rate2-needed_rate2+(Optimal_rate1-needed_rate1).*(needed_rate1<=Optimal_rate2) ).*...
%             (needed_rate2<=Optimal_rate3) ).*(needed_rate3<=Optimal_rate4) ).*(Optimal_rate4<=Perfect_rate4))...
%             ./(N+sum(Second_T)+sum(Third_T)+sum(Fourth_T));
%
%         throughput_harq_temp(4,ind_delta)=sum( Optimal_rate1.*(Optimal_rate1<=Perfect_rate1)+ (Optimal_rate1>Perfect_rate1).*...
%             Optimal_rate1.*(Optimal_rate1<=Perfect_rate1+Perfect_rate2) +(Optimal_rate1>Perfect_rate1)...
%             .*(Optimal_rate1>Perfect_rate1+Perfect_rate2).*Optimal_rate1.*(Optimal_rate1<=Perfect_rate1+Perfect_rate2+Perfect_rate3)...
%             +(Optimal_rate1>Perfect_rate1).*(Optimal_rate1>Perfect_rate1+Perfect_rate2).* (Optimal_rate1>Perfect_rate1+Perfect_rate2+Perfect_rate3).*...
%             Optimal_rate1.*(Optimal_rate1<=Perfect_rate1+Perfect_rate2+Perfect_rate3+Perfect_rate4) )...
%             ./( N+sum(Second_T)+sum(Third_T_C)+sum(Fourth_T_C) );
%         end
%     end
%
%     for ind_k=1:K
%         [throughput_lharq(ind_k,ind_snr),ind] = max(throughput_lharq_temp(ind_k,:));
%         delta_lharq(ind_k,ind_snr)=delta_set(ind);
%         outage(ind_k,ind_snr)=outage_temp(ind);
%         [throughput_harq(ind_k,ind_snr),ind] = max(throughput_harq_temp(ind_k,:));
%         delta_harq(ind_k,ind_snr)=delta_set(ind);
%     end


% end
%save LHARQ_translated_AMC_policy_16QAM_Perfect_Decoding_continuous_R_dop_0.075rho_0.9_V_leger.mat  dop rho SNRdB ...
%   throughput_lharq throughput_harq outage snr_range policy delta_harq delta_lharq

