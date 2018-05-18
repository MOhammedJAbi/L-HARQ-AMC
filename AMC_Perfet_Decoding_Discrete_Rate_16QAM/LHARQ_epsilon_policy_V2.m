%clc;
clear;

SNRdB=[-5:30];

dop=0.05
rho=(besselj(0,2*pi*dop)).^2

% pdf joint, av_snr=channel average SNR
pdf_snr=@(x,av_snr) exp(-x./av_snr)./av_snr;
% SNR cdf
cdf_snr=@(x,av_snr) 1-exp(-x./av_snr);
% optimal threshold

load MI_16QAM.mat
MI_16QAM=[0 MI_16QAM(1:end-1)' 4];
th=[0 10.^((SNR(1:end-1))'./10) 1e10];


% plot(th(1:end-1),MI_16QAM(1:end-1))

Qam_th=@(R) interp1(MI_16QAM,th,R)
Qam_r=@(gamma) interp1(th,MI_16QAM,gamma)

f_threshold1=@(x,R,av_snr) 1-marcumq(sqrt(rho*2.*x./((1-rho).*av_snr)),sqrt(2.*Qam_th(R)./((1-rho).*av_snr)));

% Legend integral
lint =2e2;
[x0, w] = GaussLegendre(lint);

epsilon_set=logspace(-1.7,-0.01,30);
N=1e7;
K=4;



%% Computation of throughput
for ind_snr=28:length(SNRdB)
    ind_snr
    snr=10.^(SNRdB(ind_snr)./10);
    max_rate=4;
    Dx=snr.*8;
    x= (x0+1)*Dx/2;
    
    y=([linspace(1e-5,MI_16QAM(end-11),50) linspace(MI_16QAM(end-10),4-1e-5,3e2) ] )';
    
    [xmat,ymat] = meshgrid(x,y);
    f=f_threshold1(xmat,ymat,snr);
    
    figure(1)
    j=10
    x(j)
    semilogy(y,f(:,j) )
    
     figure(2)
    j=80
    x(j)
    semilogy(y,f(:,j) )
    
    sigma=sqrt((1-rho).*snr./2);
    eR_m = normrnd(0,sigma,K,N);
    eI_m = normrnd(0,sigma,K,N);
    heR_m = normrnd(0,sqrt(snr./2),K,N);
    heI_m = normrnd(0,sqrt(snr./2),K,N);
    
    for ind_epsilon=1:length(epsilon_set)
        
        for j=1:length(x)
            ff=exp(f(:,j));
            [v,ind]=unique(ff,'last');
            ee=exp(epsilon_set(ind_epsilon));
            policy_tem(j)=interp1(v,y(ind),ee);
            %             [v,ind]=unique(f(:,j),'last');
            %
            %             policy_tem(j)=interp1(v,y(ind),epsilon_set(ind_epsilon));
        end
        
        
        policy_m=policy_tem;
        policy_m=[policy_m(1) policy_m policy_m(end)];
        snr_m=x';
        snr_m=[0 snr_m 1e100];
        %         epsilon_set(ind_epsilon)
        %         figure(1)
        %         semilogy(x,f_threshold1(x,policy_tem',snr) )
        %
        %                 figure(2)
        %                  semilogx(x',policy_tem)
        
        %
        %% K=1
        
        tt=1;
        eR = eR_m(tt,:);
        eI = eI_m(tt,:);
        heR = heR_m(tt,:);
        heI =heI_m(tt,:);
        gammae=heR.^2+ heI.^2;
        gamma=(sqrt(rho).*heR-eR).^2+(sqrt(rho).*heI-eI).^2;
        Optimal_rate1=interp1(snr_m,policy_m,gammae);
        Perfect_rate1=Qam_r(gamma);%log2(1+gamma);
        needed_rate1=max(0,Optimal_rate1-Perfect_rate1);
        
        throughput_lharq_temp(1,ind_epsilon)=sum(Optimal_rate1.*(Optimal_rate1<=Perfect_rate1))./N;
        throughput_harq_temp(1,ind_epsilon)=throughput_lharq_temp(1,ind_epsilon);
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
        Perfect_rate2=Qam_r(gamma2); %log2(1+gamma2);
        needed_rate2=max(0,Optimal_rate2-Perfect_rate2);
        Second_T=(Optimal_rate1>Perfect_rate1).*(needed_rate1<=Optimal_rate2);
        throughput_lharq_temp(2,ind_epsilon)=sum( Optimal_rate1.*(Optimal_rate1<=Perfect_rate1) + (Optimal_rate1>Perfect_rate1).*(needed_rate1<=Optimal_rate2).*...
            ( Optimal_rate2+(Optimal_rate1-needed_rate1).*(needed_rate1<=Optimal_rate2) ) .* ( Optimal_rate2<=Perfect_rate2) )...
            ./( N+sum( Second_T) );
        throughput_harq_temp(2,ind_epsilon)=sum( Optimal_rate1.*(Optimal_rate1<=Perfect_rate1) + (Optimal_rate1>Perfect_rate1).*...
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
        Perfect_rate3=Qam_r(gamma3);%log2(1+gamma3);
        needed_rate3=max(0,Optimal_rate3-Perfect_rate3);
        Third_T=(Optimal_rate1>Perfect_rate1).*(Optimal_rate2>Perfect_rate2).*(needed_rate1<=Optimal_rate2).*(needed_rate2<=Optimal_rate3);
        Third_T_C=(Optimal_rate1>Perfect_rate1).*(Optimal_rate1>Perfect_rate1+Perfect_rate2);
        
        throughput_lharq_temp(3,ind_epsilon)=sum( Optimal_rate1.*(Optimal_rate1<=Perfect_rate1) + (Optimal_rate1>Perfect_rate1).*(needed_rate1<=Optimal_rate2).*...
            ( Optimal_rate2+(Optimal_rate1-needed_rate1).*(needed_rate1<=Optimal_rate2) ) .* ( Optimal_rate2<=Perfect_rate2) +...
            (Optimal_rate1>Perfect_rate1).*(Optimal_rate2>Perfect_rate2).*(needed_rate1<=Optimal_rate2).*(needed_rate2<=Optimal_rate3).*...
            ( Optimal_rate3+ ( Optimal_rate2-needed_rate2+(Optimal_rate1-needed_rate1).*(needed_rate1<=Optimal_rate2) ).*...
            (needed_rate2<=Optimal_rate3) ).*(Optimal_rate3<=Perfect_rate3) )  ...
            ./(N+sum(Second_T)+sum(Third_T));
        
        throughput_harq_temp(3,ind_epsilon)=sum( Optimal_rate1.*(Optimal_rate1<=Perfect_rate1)+ (Optimal_rate1>Perfect_rate1).*...
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
        Perfect_rate4=Qam_r(gamma4) ;%log2(1+gamma4);
        needed_rate4=max(0,Optimal_rate4-Perfect_rate4);
        Fourth_T=(Optimal_rate1>Perfect_rate1).*(Optimal_rate2>Perfect_rate2).*(Optimal_rate3>Perfect_rate3).*(needed_rate1<=Optimal_rate2)...
            .*(needed_rate2<=Optimal_rate3).*(needed_rate3<=Optimal_rate4);
        Fourth_T_C=(Optimal_rate1>Perfect_rate1).*(Optimal_rate1>Perfect_rate1+Perfect_rate2).*(Optimal_rate1>Perfect_rate1+Perfect_rate2+Perfect_rate3);
        
        throughput_lharq_temp(4,ind_epsilon)=sum( Optimal_rate1.*(Optimal_rate1<=Perfect_rate1) + (Optimal_rate1>Perfect_rate1).*(needed_rate1<=Optimal_rate2).*...
            ( Optimal_rate2+(Optimal_rate1-needed_rate1).*(needed_rate1<=Optimal_rate2) ) .* ( Optimal_rate2<=Perfect_rate2) +...
            (Optimal_rate1>Perfect_rate1).*(Optimal_rate2>Perfect_rate2).*(needed_rate1<=Optimal_rate2).*(needed_rate2<=Optimal_rate3).*...
            ( Optimal_rate3+ ( Optimal_rate2-needed_rate2+(Optimal_rate1-needed_rate1).*(needed_rate1<=Optimal_rate2) ).*...
            (needed_rate2<=Optimal_rate3) ).*(Optimal_rate3<=Perfect_rate3)+  ...
            (Optimal_rate1>Perfect_rate1).*(Optimal_rate2>Perfect_rate2).*(Optimal_rate3>Perfect_rate3).*(needed_rate1<=Optimal_rate2).*(needed_rate2<=Optimal_rate3).*(needed_rate3<=Optimal_rate4).*...
            ( Optimal_rate4+ ( Optimal_rate3-needed_rate3+( Optimal_rate2-needed_rate2+(Optimal_rate1-needed_rate1).*(needed_rate1<=Optimal_rate2) ).*...
            (needed_rate2<=Optimal_rate3) ).*(needed_rate3<=Optimal_rate4) ).*(Optimal_rate4<=Perfect_rate4))...
            ./(N+sum(Second_T)+sum(Third_T)+sum(Fourth_T));
        
        throughput_harq_temp(4,ind_epsilon)=sum( Optimal_rate1.*(Optimal_rate1<=Perfect_rate1)+ (Optimal_rate1>Perfect_rate1).*...
            Optimal_rate1.*(Optimal_rate1<=Perfect_rate1+Perfect_rate2) +(Optimal_rate1>Perfect_rate1)...
            .*(Optimal_rate1>Perfect_rate1+Perfect_rate2).*Optimal_rate1.*(Optimal_rate1<=Perfect_rate1+Perfect_rate2+Perfect_rate3)...
            +(Optimal_rate1>Perfect_rate1).*(Optimal_rate1>Perfect_rate1+Perfect_rate2).* (Optimal_rate1>Perfect_rate1+Perfect_rate2+Perfect_rate3).*...
            Optimal_rate1.*(Optimal_rate1<=Perfect_rate1+Perfect_rate2+Perfect_rate3+Perfect_rate4) )...
            ./( N+sum(Second_T)+sum(Third_T_C)+sum(Fourth_T_C) );
    end
    snr_range(ind_snr,:)=x;
    %%
    ii=1;
    [throughput_lharq(ii,ind_snr),I]=max(throughput_lharq_temp(ii,:));
    optimal_epsilon_lharq(ii,ind_snr)=epsilon_set(I);

    for j=1:length(x)
        [v,ind]=unique(f(:,j),'last');
        policy_tem(j)=interp1(v,y(ind),epsilon_set(I));
    end
    policy_lharq(ii,ind_snr,:)=policy_tem;
    
    policy_m=policy_tem;
    policy_m=[policy_m(1) policy_m policy_m(end)];
    snr_m=x';
    snr_m=[0 snr_m 1e100];
    tt=1;
    eR = eR_m(tt,:);
    eI = eI_m(tt,:);
    heR = heR_m(tt,:);
    heI =heI_m(tt,:);
    gammae=heR.^2+ heI.^2;
    gamma=(sqrt(rho).*heR-eR).^2+(sqrt(rho).*heI-eI).^2;
    Optimal_rate1=interp1(snr_m,policy_m,gammae);
    Perfect_rate1=Qam_r(gamma);%log2(1+gamma);
    needed_rate1=max(0,Optimal_rate1-Perfect_rate1);
    optimal_epsilon_lharq_p(ii,ind_snr)=sum((Optimal_rate1>Perfect_rate1))./N;
    
%     figure(3)
%     plot(x,squeeze(policy_lharq(ii,ind_snr,:)))
%     figure(4)
%     semilogy(x,f_threshold1(x,squeeze(policy_lharq(ii,ind_snr,:)),snr) )
    %        1
    
    
    %
    ii=2;
    [throughput_lharq(ii,ind_snr),I]=max(throughput_lharq_temp(ii,:));
    optimal_epsilon_lharq(ii,ind_snr)=epsilon_set(I);
    for j=1:length(x)
        [v,ind]=unique(f(:,j),'last');
        policy_tem(j)=interp1(v,y(ind),epsilon_set(I));
    end
    policy_lharq(ii,ind_snr,:)=policy_tem;
    
    [throughput_harq(ii,ind_snr),I]=max(throughput_harq_temp(ii,:));
    optimal_epsilon_harq(ii,ind_snr)=epsilon_set(I);
     policy_m=policy_tem;
    policy_m=[policy_m(1) policy_m policy_m(end)];
    snr_m=x';
    snr_m=[0 snr_m 1e100];
    tt=1;
    eR = eR_m(tt,:);
    eI = eI_m(tt,:);
    heR = heR_m(tt,:);
    heI =heI_m(tt,:);
    gammae=heR.^2+ heI.^2;
    gamma=(sqrt(rho).*heR-eR).^2+(sqrt(rho).*heI-eI).^2;
    Optimal_rate1=interp1(snr_m,policy_m,gammae);
    Perfect_rate1=Qam_r(gamma);%log2(1+gamma);
    needed_rate1=max(0,Optimal_rate1-Perfect_rate1);
    optimal_epsilon_lharq_p(ii,ind_snr)=sum((Optimal_rate1>Perfect_rate1))./N;
    
    for j=1:length(x)
        [v,ind]=unique(f(:,j),'last');
        policy_tem(j)=interp1(v,y(ind),epsilon_set(I));
    end
    policy_harq(ii,ind_snr,:)=policy_tem;
    %%
    ii=3;
    [throughput_lharq(ii,ind_snr),I]=max(throughput_lharq_temp(ii,:));
    optimal_epsilon_lharq(ii,ind_snr)=epsilon_set(I);
    for j=1:length(x)
        [v,ind]=unique(f(:,j),'last');
        policy_tem(j)=interp1(v,y(ind),epsilon_set(I));
    end
    policy_lharq(ii,ind_snr,:)=policy_tem;
     policy_m=policy_tem;
    policy_m=[policy_m(1) policy_m policy_m(end)];
    snr_m=x';
    snr_m=[0 snr_m 1e100];
    tt=1;
    eR = eR_m(tt,:);
    eI = eI_m(tt,:);
    heR = heR_m(tt,:);
    heI =heI_m(tt,:);
    gammae=heR.^2+ heI.^2;
    gamma=(sqrt(rho).*heR-eR).^2+(sqrt(rho).*heI-eI).^2;
    Optimal_rate1=interp1(snr_m,policy_m,gammae);
    Perfect_rate1=Qam_r(gamma);%log2(1+gamma);
    needed_rate1=max(0,Optimal_rate1-Perfect_rate1);
    optimal_epsilon_lharq_p(ii,ind_snr)=sum((Optimal_rate1>Perfect_rate1))./N;
    
    
    [throughput_harq(ii,ind_snr),I]=max(throughput_harq_temp(ii,:));
    optimal_epsilon_harq(ii,ind_snr)=epsilon_set(I);
    for j=1:length(x)
        [v,ind]=unique(f(:,j),'last');
        policy_tem(j)=interp1(v,y(ind),epsilon_set(I));
    end
    policy_harq(ii,ind_snr,:)=policy_tem;
    %%
    ii=4;
    [throughput_lharq(ii,ind_snr),I]=max(throughput_lharq_temp(ii,:));
    optimal_epsilon_lharq(ii,ind_snr)=epsilon_set(I);
    for j=1:length(x)
        [v,ind]=unique(f(:,j),'last');
        policy_tem(j)=interp1(v,y(ind),epsilon_set(I));
    end
    policy_lharq(ii,ind_snr,:)=policy_tem;
     policy_m=policy_tem;
    policy_m=[policy_m(1) policy_m policy_m(end)];
    snr_m=x';
    snr_m=[0 snr_m 1e100];
    tt=1;
    eR = eR_m(tt,:);
    eI = eI_m(tt,:);
    heR = heR_m(tt,:);
    heI =heI_m(tt,:);
    gammae=heR.^2+ heI.^2;
    gamma=(sqrt(rho).*heR-eR).^2+(sqrt(rho).*heI-eI).^2;
    Optimal_rate1=interp1(snr_m,policy_m,gammae);
    Perfect_rate1=Qam_r(gamma);%log2(1+gamma);
    needed_rate1=max(0,Optimal_rate1-Perfect_rate1);
    optimal_epsilon_lharq_p(ii,ind_snr)=sum((Optimal_rate1>Perfect_rate1))./N;
    
    
    [throughput_harq(ii,ind_snr),I]=max(throughput_harq_temp(ii,:));
    optimal_epsilon_harq(ii,ind_snr)=epsilon_set(I);
    for j=1:length(x)
        [v,ind]=unique(f(:,j),'last');
        policy_tem(j)=interp1(v,y(ind),epsilon_set(I));
    end
    policy_harq(ii,ind_snr,:)=policy_tem;
    
    
    
    
end

save LHARQ_outage_policy_16QAM_Perfect_Decoding_continuous_R_dop_0.05_rho_0.95_V_leger_V2.mat dop rho SNRdB snr_range ...
  throughput_lharq optimal_epsilon_lharq policy_lharq throughput_harq  optimal_epsilon_harq policy_harq optimal_epsilon_lharq_p



%save(['LHARQ_outage_policy_16QAM_Perfect_Decoding_continuous_R_dop_',num2str(dop),'_rho_',num2str(rho),'.mat'])