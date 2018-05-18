clc;
clear;
%clf;


SNRdB=[-5:30];
% correlation coefficient
dop=0.1;
rho=(besselj(0,2*pi*dop)).^2;
erg=@(x,av_snr) log2(1+x).*exp(-x./av_snr)./av_snr;

lambdaMin=1e-2;


K=2;

%% Computation of throughput
for ind_snr=1:length(SNRdB)
    ind_snr
    snr=10.^(SNRdB(ind_snr)./10);
    ergodic_capacity(ind_snr)=integral(@(x) erg(x,snr),0,inf);
    max_rate=max(2.*ergodic_capacity(ind_snr),15);
    lambdaMax=max(1.2.*ergodic_capacity(ind_snr),14);
    lambda_set=linspace(lambdaMin,lambdaMax,30);
    lint = 1e2;
    [x0, w] = GaussLegendre(lint);
    Dx=15.*snr;
    x = (x0+1)*Dx/2;
    length_rate=50;
    
    need_set=linspace(0,max_rate,length_rate);
    rate_set=linspace(0.1,max_rate,length_rate);
    [Mat_rate,Mat_g]=meshgrid(rate_set,x);
    [Mat_need,Mat_rate2]=meshgrid(need_set,rate_set');
    snr_m=x';
    snr_m=[0 snr_m 1e100];
    [Mat_rate3,Mat_snr]=meshgrid(rate_set,snr_m');
    
    N=1e6;
    sigma=sqrt((1-rho).*snr./2);
    eR_m = normrnd(0,sigma,K,N);
    eI_m = normrnd(0,sigma,K,N);
    heR_m = normrnd(0,sqrt(snr./2),K,N);
    heI_m = normrnd(0,sqrt(snr./2),K,N);
    %     tt=1;
    %     eR = eR_m(tt,:);
    %     eI = eI_m(tt,:);
    %     heR = heR_m(tt,:);
    %     heI =heI_m(tt,:);
    %     gammae=heR.^2+ heI.^2;
    %     gamma=(sqrt(rho).*heR-eR).^2+(sqrt(rho).*heI-eI).^2;
    %     gamma=log2(1+gamma);
    %     moy=sum(gamma)./N;
    
    pdf_snr=@(g) exp(-g./snr)./snr;
    pdf_cond=@(g,ge) exp(-(g+rho.*ge-2*sqrt(rho.*g.*ge))./((1-rho).*snr)).*besseli(0,2.*sqrt(rho.*g.*ge)./((1-rho).*snr),1)./((1-rho).*snr);
    
    
    for ind_lambda=1:length(lambda_set)
        
        lambda=lambda_set(ind_lambda);
        V_2=@(ge2,R2,R1) (R2+R1).*marcumq(sqrt(rho*2.*ge2./((1-rho).*snr)),sqrt(2.*(2.^R2-1)./((1-rho).*snr)));
        V_2c=@(ge2,R2,R1,n) (R2-n+R1).*(R2>=n).*marcumq(sqrt(rho*2.*ge2./((1-rho).*snr)),sqrt(2.*(2.^R2-1)./((1-rho).*snr)))-lambda.*(R2>=n);
        MParity_Rate2=zeros(length(need_set),length(x)+2);
        
        tic
        for i=1:length(rate_set)
            [val0,ind]=max(V_2(Mat_g,Mat_rate,rate_set(i))');
            
            Parity_Rate2(i,:)=rate_set(ind);
            for j=1:length(need_set)
                val=V_2c(x,rate_set(ind)',rate_set(i),need_set(j));
                xtemp=val.*pdf_snr(x);
                V_av2(j,i)=sum(xtemp.*w)*Dx/2;
            end
        end
        toc
        
        
        MParity_Rate2=[Parity_Rate2(:,1) Parity_Rate2 Parity_Rate2(:,end)];
        FInterpolation = scatteredInterpolant(Mat_need(:),Mat_rate2(:),V_av2(:));
        policyInterpolation = scatteredInterpolant(Mat_rate3(:),Mat_snr(:),MParity_Rate2(:));
        
        J1=@(R1,ge1,g1) FInterpolation(R1-log2(1+g1),R1).*pdf_cond(g1,ge1);
        V_1=@(R1,ge1) R1.*marcumq(sqrt(rho*2.*ge1./((1-rho).*snr)),sqrt(2.*(2.^R1-1)./((1-rho).*snr)));%+integral(@(g1) FInterpolation(R-log2(1+g1),log2(1+g1)).*pdf_cond(g1,ge1),0, 2.^R-1 );
        
        tic
        for i=1:length(x)
            for j=1:length(rate_set)
                th_d=2.^rate_set(j)-1;
                y = (x0+1)*th_d/2;
                val_temp=J1(rate_set(j).*ones(lint,1),x(i),y);
                temp(j)=sum(val_temp.*w)*th_d/2+V_1(rate_set(j),x(i));
            end
            [val,ind0]=max(temp);
            V{1}(i)=val;
            Parity_Rate1(i)=rate_set(ind0);
        end
        toc
        
        policy_m=Parity_Rate1;
        policy_m=[policy_m(1) policy_m policy_m(end)];
        
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
        
        throughput_lharq_temp(1,ind_lambda)=sum(Optimal_rate1.*(Optimal_rate1<=Perfect_rate1))./N;
        throughput_harq_temp(1,ind_lambda)=throughput_lharq_temp(1,ind_lambda);
        outage(ind_lambda)=sum((Optimal_rate1>Perfect_rate1))./N;
        
        %% K=2
        
        tt=2;
        eR = eR_m(tt,:);
        eI = eI_m(tt,:);
        heR = heR_m(tt,:);
        heI =heI_m(tt,:);
        gammae2=heR.^2+ heI.^2;
        gamma2=(sqrt(rho).*heR-eR).^2+(sqrt(rho).*heI-eI).^2;
        Optimal_rate2=policyInterpolation(Optimal_rate1,gammae2);% interp1(snr_m,policy_m,gammae2);
        Perfect_rate2=log2(1+gamma2);
        needed_rate2=max(0,Optimal_rate2-Perfect_rate2);
        Second_T=(Optimal_rate1>Perfect_rate1).*(needed_rate1<=Optimal_rate2);
        throughput_lharq_temp(2,ind_lambda)=sum( Optimal_rate1.*(Optimal_rate1<=Perfect_rate1) + (Optimal_rate1>Perfect_rate1).*(needed_rate1<=Optimal_rate2).*...
            ( Optimal_rate2+(Optimal_rate1-needed_rate1).*(needed_rate1<=Optimal_rate2) ) .* ( Optimal_rate2<=Perfect_rate2) )...
            ./( N+sum( Second_T) );
        throughput_harq_temp(2,ind_lambda)=sum( Optimal_rate1.*(Optimal_rate1<=Perfect_rate1) + (Optimal_rate1>Perfect_rate1).*...
            Optimal_rate1.*(Optimal_rate1<=Perfect_rate1+Perfect_rate2) )...
            ./( N+sum( Optimal_rate1>Perfect_rate1 ) );
        
    end
    
    %%
    ii=2;
    [throughput_lharq(ii,ind_snr),I]=max(throughput_lharq_temp(ii,:));
    optimal_lambda_lharq(ii,ind_snr)=lambda_set(I);
    
    
    
    lambda=lambda_set(I);
    
    pdf_snr=@(g) exp(-g./snr)./snr;
    pdf_cond=@(g,ge) exp(-(g+rho.*ge-2*sqrt(rho.*g.*ge))./((1-rho).*snr)).*besseli(0,2.*sqrt(rho.*g.*ge)./((1-rho).*snr),1)./((1-rho).*snr);
    
    V_2=@(ge2,R2,R1) (R2+R1).*marcumq(sqrt(rho*2.*ge2./((1-rho).*snr)),sqrt(2.*(2.^R2-1)./((1-rho).*snr)));
    V_2c=@(ge2,R2,R1,n) (R2-n+R1).*(R2>=n).*marcumq(sqrt(rho*2.*ge2./((1-rho).*snr)),sqrt(2.*(2.^R2-1)./((1-rho).*snr)))-lambda.*(R2>=n);
    MParity_Rate2=zeros(length(need_set),length(x)+2);
    
    for i=1:length(rate_set)
        [val0,ind]=max(V_2(Mat_g,Mat_rate,rate_set(i))');
        
        Parity_Rate2(i,:)=rate_set(ind);
        for j=1:length(need_set)
            val=V_2c(x,rate_set(ind)',rate_set(i),need_set(j));
            xtemp=val.*pdf_snr(x);
            V_av2(j,i)=sum(xtemp.*w)*Dx/2;
        end
    end
    
    policy_lharq{2}(ind_snr,:,:)=Parity_Rate2;
    
    MParity_Rate2=[Parity_Rate2(:,1) Parity_Rate2 Parity_Rate2(:,end)];
    FInterpolation = scatteredInterpolant(Mat_need(:),Mat_rate2(:),V_av2(:));
    policyInterpolation = scatteredInterpolant(Mat_rate3(:),Mat_snr(:),MParity_Rate2(:));
    
    J1=@(R1,ge1,g1) FInterpolation(R1-log2(1+g1),R1).*pdf_cond(g1,ge1);
    V_1=@(R1,ge1) R1.*marcumq(sqrt(rho*2.*ge1./((1-rho).*snr)),sqrt(2.*(2.^R1-1)./((1-rho).*snr)));%+integral(@(g1) FInterpolation(R-log2(1+g1),log2(1+g1)).*pdf_cond(g1,ge1),0, 2.^R-1 );
    
    
    for i=1:length(x)
        for j=1:length(rate_set)
            th_d=2.^rate_set(j)-1;
            y = (x0+1)*th_d/2;
            val_temp=J1(rate_set(j).*ones(lint,1),x(i),y);
            temp(j)=sum(val_temp.*w)*th_d/2+V_1(rate_set(j),x(i));
        end
        [val,ind0]=max(temp);
        V{1}(i)=val;
        Parity_Rate1(i)=rate_set(ind0);
    end
    
    
    policy_lharq{1}(ind_snr,:)=Parity_Rate1';
    
    
    
    
end

save LHARQ_AMC_DP_policy_Perfect_Decoding_continuous_R_dop_0.1_rho_0.8_V3.mat SNRdB lambda_set optimal_lambda_lharq rho policy_lharq throughput_lharq
