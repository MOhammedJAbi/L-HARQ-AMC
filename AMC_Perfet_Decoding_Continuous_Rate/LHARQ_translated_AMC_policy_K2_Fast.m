clear;
%clc;
%clf;



[FInterpolation1,FInterpolation2,FInterpolation3,FInterpolation4] = WEP_space(1);

epsilon_Set=[logspace(-2,-1.1,10)  linspace(0.08,0.95,10)];

%% determination of back track rate policies
[parity_policy,gamma_interval] = parity_policy(epsilon_Set);

load AMC_policy_16QAM_Turbo_code_dop_0.05_rho_0.95.mat

SNRdB=[-5:30];
K=2;
N_realization=1e5;

for ind_snr=1:1:length(SNRdB)
    
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
            
    delta_set=logspace(0,log10(snr_range(ind_snr,end-1)),20);
    
    for ind_delta=1:length(delta_set)
       [ind_snr ind_delta]
        snr_mm=snr_m./delta_set(ind_delta);
        
        for ind_epsilon=1:length(epsilon_Set)
            [ind_snr ind_epsilon]
            counter=0;
            reward=0;
            number_ack_1=0;
            while counter<N_realization-1
                counter=counter+1;
                rate_round_1=interp1(snr_mm,policy_m,gammae(counter),'nearest');
                wep_round_1=WEP_function_Fast(rate_round_1,0,gamma(counter),FInterpolation1,FInterpolation2,FInterpolation3,FInterpolation4);
                ack_round_1=binornd(1, 1-wep_round_1);
                    if ack_round_1==1
                        reward=reward+rate_round_1;
                        number_ack_1=number_ack_1+1;
                    else
                        index_rate_round_1 = index_rate(rate_round_1);
                        BK_rate_round_1 = interp1(gamma_interval(index_rate_round_1,:)...
                            ,squeeze(parity_policy(index_rate_round_1,ind_epsilon,:)),gamma(counter),'previous');
                        rate_round_2=interp1(snr_mm,policy_m,gammae(counter+1),'nearest');
                        if BK_rate_round_1<=rate_round_2
                            BK_wep_round_1=WEP_function_Fast(rate_round_1,BK_rate_round_1,gamma(counter),FInterpolation1,FInterpolation2,FInterpolation3,FInterpolation4);
                            BK_wep_round_1=min(max(BK_wep_round_1/wep_round_1,0),1);
                            BK_ack_round_1=binornd(1, 1-BK_wep_round_1);
                            wep_round_2=WEP_function_Fast(rate_round_2,0,gamma(counter+1),FInterpolation1,FInterpolation2,FInterpolation3,FInterpolation4);
                            ack_round_2=binornd(1, 1-wep_round_2);
                            reward=reward+(rate_round_2+(rate_round_1-BK_rate_round_1)*BK_ack_round_1)...
                                *ack_round_2;
                            counter=counter+1;
                        end
                    end
            end
            throughput_lharq_epsilon(ind_epsilon,ind_delta,ind_snr)= reward/counter;
        end
        [throughput_lharq_delta(ind_delta,ind_snr),ii]=max(throughput_lharq_epsilon(:,ind_delta,ind_snr));
        optimal_epsilon_delta(ind_delta,ind_snr)=epsilon_Set(ii(1));
    end
    [throughput_lharq(ind_snr),jj]=max(throughput_lharq_delta(:,ind_snr));
    optimal_delta(ind_snr)=delta_set(jj);
    optimal_epsilon(ind_snr)=optimal_epsilon_delta(jj,ind_snr);
end


save LHARQ_translated_AMC_policy_16QAM_turbo_code_K_2_dop_0.05_rho_0.95.mat  dop rho SNRdB throughput_lharq_delta optimal_epsilon_delta...
  throughput_lharq optimal_delta optimal_epsilon



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

