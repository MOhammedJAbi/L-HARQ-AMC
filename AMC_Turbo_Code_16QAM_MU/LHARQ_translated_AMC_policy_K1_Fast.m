clear;clc;clf;

[FInterpolation1,FInterpolation2,FInterpolation3,FInterpolation4] = WEP_space(1);

load AMC_policy_16QAM_Turbo_code_dop_0.05rho_0.95156.mat
throughput_amc=0;
L_set=[10:10:100];
path_loss=3;
SNR_function=@(d) d.^(-path_loss);
K=1;
N_realization=1e4;
N=1e3;
for ind_L=1:length(L_set)
    ind_L
    L=L_set(ind_L);
    accumulated_reward=1e-6.*ones(1,L);
    counter=1;
    while counter<N_realization
        d=sqrt(unifrnd(0,1,1,L));
        d(d<0.1)=0.1;
        snr=SNR_function(d);
        sigma=sqrt((1-rho).*snr./2);
        gammae=exprnd(snr);
        phi = unifrnd(0,2.*pi,N,L);
        eR = normrnd(zeros(N,L),repmat(sigma,N,1),N,L);
        eI = normrnd(zeros(N,L),repmat(sigma,N,1),N,L);
        gamma=rho.*repmat(gammae,N,1)+eI.^2+eR.^2-2.*sqrt(rho).*...
            eR.*sqrt(repmat(gammae,N,1)).*cos(phi)-2.*sqrt(rho).*eI.*sqrt(repmat(gammae,N,1)).*sin(phi);
        %
        ind_SNRdB=interp1(SNRdB,[1:1:36],10.*log10(snr),'nearest');
        expected_throughput=0;
        for user=1:1:L
            rate_round_1=interp1(snr_range(ind_SNRdB(user),:),policy(ind_SNRdB(user),:),gammae(user),'nearest');
            wep_round_1=WEP_function_Fast_vector(rate_round_1,0,gamma(:,user),FInterpolation1,FInterpolation2,FInterpolation3,FInterpolation4);
            ack_round_1=binornd(1, 1-wep_round_1);
            expected_throughput(user)=rate_round_1.*sum(ack_round_1);
        end
        [Y,J]=max(expected_throughput./accumulated_reward);
        phi_one = unifrnd(0,2.*pi,1,1);
        eR_one = normrnd(0,sigma(J),1,1);
        eI_one = normrnd(0,sigma(J),1,1);
        gamma_one = rho.*gammae(J)+eI_one.^2+eR_one.^2-2.*sqrt(rho).*...
            eR_one.*sqrt(gammae(J)).*cos(phi_one)-2.*sqrt(rho).*eI_one.*sqrt(gammae(J)).*sin(phi_one);
        rate_round_1_one = interp1(snr_range(ind_SNRdB(J),:),policy(ind_SNRdB(J),:),gammae(J),'nearest');
        wep_round_1_one = WEP_function_Fast_vector(rate_round_1_one,0,gamma_one,FInterpolation1,FInterpolation2,FInterpolation3,FInterpolation4);
        ack_round_1_one=binornd(1, 1-wep_round_1_one);
        accumulated_reward(J)=accumulated_reward(J)+rate_round_1_one*ack_round_1_one;
        counter=counter+1;
    end
    throughput_amc(ind_L)= sum(accumulated_reward-1e-6.*ones(1,L))/counter;
end


save MU_LHARQ_translated_AMC_policy_16QAM_turbo_code_K_1_dop_0.05_rho_0.95.mat  dop rho SNRdB throughput_amc L_set



