clc;
clear;
%clf;


SNRdB=[-5:30];
% correlation coefficient
dop=0.1
rho=(besselj(0,2*pi*dop)).^2
erg=@(x,av_snr) log2(1+x).*exp(-x./av_snr)./av_snr;



lambdaMin=1e-1;
lambdaMax=14;
lambda_set=linspace(lambdaMin,lambdaMax,30)

K=2;

%% Computation of throughput
for ind_snr=15   %:length(SNRdB)
    ind_snr
    snr=10.^(SNRdB(ind_snr)./10);
    ergodic_capacity(ind_snr)=integral(@(x) erg(x,snr),0,inf);
    max_rate=max(2.*ergodic_capacity(ind_snr),15);
    lint = 55;
    [x0, w] = GaussLegendre(lint);
    Dx=15.*snr;
    x = (x0+1)*Dx/2;
    length_rate=50;
    
    need_set=linspace(0,max_rate,length_rate);
    [Mat_need,Mat_g]=meshgrid(need_set,x);
    
    rate_set=linspace(0.1,max_rate,length_rate);
    [Mat_rate,Mat_g]=meshgrid(rate_set,x);
    
    N=1e6;
    snr_m=x';
    snr_m=[0 snr_m 1e100];
    [Mat_R2,Mat_snr]=meshgrid(need_set,snr_m');
    sigma=sqrt((1-rho).*snr./2);
    eR_m = normrnd(0,sigma,K,N);
    eI_m = normrnd(0,sigma,K,N);
    heR_m = normrnd(0,sqrt(snr./2),K,N);
    heI_m = normrnd(0,sqrt(snr./2),K,N);
    tt=1;
    eR = eR_m(tt,:);
    eI = eI_m(tt,:);
    heR = heR_m(tt,:);
    heI =heI_m(tt,:);
    gammae=heR.^2+ heI.^2;
    gamma=(sqrt(rho).*heR-eR).^2+(sqrt(rho).*heI-eI).^2;
    gamma=log2(1+gamma);
    moy=sum(gamma)./N;
    
    pdf_snr=@(g) exp(-g./snr)./snr;
    pdf_cond=@(g,ge) exp(-(g+rho.*ge-2*sqrt(rho.*g.*ge))./((1-rho).*snr)).*besseli(0,2.*sqrt(rho.*g.*ge)./((1-rho).*snr),1)./((1-rho).*snr);
    
    
    for ind_lambda=3:length(lambda_set)
        lambda=lambda_set(ind_lambda);
        
        V_2=@(R,n,ge2,moy) (R+moy).*(R>=n).*marcumq(sqrt(rho*2.*ge2./((1-rho).*snr)),sqrt(2.*(2.^R-1)./((1-rho).*snr))).*(n>0)-lambda.*(n>0).*(R>=n);
        V_2c=@(R,n,ge2,g1) (R+log2(1+g1)).*(R>=n).*marcumq(sqrt(rho*2.*ge2./((1-rho).*snr)),sqrt(2.*(2.^R-1)./((1-rho).*snr))).*(n>0)-lambda.*(n>0).*(R>=n);
        MParity_Rate2=zeros(length(need_set),length(x)+2);
        tic
        
        for i=1:length(need_set)
            
            %  tic
            % [val0,ind]=max((V_2(Mat_R,R_set(i),Mat_g))');
            [val0,ind]=max(V_2(Mat_rate,need_set(i),Mat_g,moy)');
            Parity_Rate2(i,:)=rate_set(ind);
            for j=1:length(x)
                val=V_2c(rate_set(ind)',need_set(i),x,x(j));
                xtemp=val.*pdf_snr(x);
                V_av2(i,j)=sum(xtemp.*w)*Dx/2;
            end
            %  toc
        end
        MParity_Rate2=[Parity_Rate2(:,1) Parity_Rate2 Parity_Rate2(:,end)];
        toc
        FInterpolation = scatteredInterpolant(Mat_need(:),Mat_g(:),V_av2(:));
        
        policyInterpolation = scatteredInterpolant(Mat_R2(:),Mat_snr(:),MParity_Rate2(:));
        
        J1=@(R,ge1,g1) FInterpolation(R-log2(1+g1),g1).*pdf_cond(g1,ge1);
        
        V_1=@(R,ge1) R.*marcumq(sqrt(rho*2.*ge1./((1-rho).*snr)),sqrt(2.*(2.^R-1)./((1-rho).*snr)));%+integral(@(g1) FInterpolation(R-log2(1+g1),log2(1+g1)).*pdf_cond(g1,ge1),0, 2.^R-1 );
        tic
        for i=1:length(x)
            for j=1:length(rate_set)
                th_d=2.^rate_set(j)-1;
                y = (x0+1)*th_d/2;
                val_temp=J1(rate_set(j),x(i),y);
                temp(j)=sum(val_temp.*w)*th_d/2+V_1(rate_set(j),x(i));
            end
            [val,ind0]=max(temp);
            V{1}(i)=val;
            Parity_Rate1(i)=rate_set(ind0);
        end
        % K=1
        
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
        Optimal_rate2=policyInterpolation(needed_rate1,gammae2);% interp1(snr_m,policy_m,gammae2);
        Perfect_rate2=log2(1+gamma2);
        needed_rate2=max(0,Optimal_rate2-Perfect_rate2);
        Second_T=(Optimal_rate1>Perfect_rate1).*(needed_rate1<=Optimal_rate2);
        throughput_lharq_temp(2,ind_lambda)=sum( Optimal_rate1.*(Optimal_rate1<=Perfect_rate1) + (Optimal_rate1>Perfect_rate1).*(needed_rate1<=Optimal_rate2).*...
            ( Optimal_rate2+(Optimal_rate1-needed_rate1).*(needed_rate1<=Optimal_rate2) ) .* ( Optimal_rate2<=Perfect_rate2) )...
            ./( N+sum( Second_T) )
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
    
    V_2=@(R,n,ge2) R.*(R>=n).*marcumq(sqrt(rho*2.*ge2./((1-rho).*snr)),sqrt(2.*(2.^R-1)./((1-rho).*snr))).*(n>0)-lambda.*(n>0).*(R>=n);
    V_2c=@(R,n,ge2,g1) (R+log2(1+g1)).*(R>=n).*marcumq(sqrt(rho*2.*ge2./((1-rho).*snr)),sqrt(2.*(2.^R-1)./((1-rho).*snr))).*(n>0)-lambda.*(n>0).*(R>=n);
    MParity_Rate2=zeros(length(R_set),length(x)+2);
    tic
    for i=1:length(R_set)
        %  tic
        %[val0,ind]=max((V_2(Mat_R,R_set(i),Mat_g))');
        [val0,ind]=max(V_2c(Mat_R,R_set(i),Mat_g,moy)');
        Parity_Rate2(i,:)=R_set(ind);
        for j=1:length(x)
            val=V_2c(R_set(ind)',R_set(i),x,x(j));
            xtemp=val.*pdf_snr(x);
            V_av2(i,j)=sum(xtemp.*w)*Dx/2;
        end
        %  toc
    end
    policy_lharq{2}(ind_snr,:,:)=Parity_Rate2;
    
    MParity_Rate2=[Parity_Rate2(:,1) Parity_Rate2 Parity_Rate2(:,end)];
    toc
    FInterpolation = scatteredInterpolant(Mat_R(:),Mat_g(:),V_av2(:));
    policyInterpolation = scatteredInterpolant(Mat_R2(:),Mat_snr(:),MParity_Rate2(:));
    
    J1=@(R,ge1,g1) FInterpolation(R-log2(1+g1),g1).*pdf_cond(g1,ge1);
    
    V_1=@(R,ge1) R.*marcumq(sqrt(rho*2.*ge1./((1-rho).*snr)),sqrt(2.*(2.^R-1)./((1-rho).*snr)));%+integral(@(g1) FInterpolation(R-log2(1+g1),log2(1+g1)).*pdf_cond(g1,ge1),0, 2.^R-1 );
    tic
    for i=1:length(x)
        for j=1:length(R_set)
            th_d=2.^R_set(j)-1;
            y = (x0+1)*th_d/2;
            val_temp=J1(R_set(j),x(i),y);
            temp(j)=sum(val_temp.*w)*th_d/2+V_1(R_set(j),x(i));
        end
        [val,ind]=max(temp);
        V{1}(i)=val;
        Parity_Rate1(i)=R_set(ind);
    end
    
    policy_lharq{1}(ind_snr,:)=Parity_Rate1';
    
    
    
    
end