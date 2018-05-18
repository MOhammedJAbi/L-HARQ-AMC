function [test_lambda]=find_lambda_DP(lambda_set,snr,max_rate,K,rho)



Dx=snr.*15;
lint = 1e2;
[x0, w] = GaussLegendre(lint);
x = (x0+1)*Dx/2;
R_set=linspace(0,max_rate,50);
[Mat_R,Mat_g]=meshgrid(R_set,x);

N=1e6;
snr_m=x';
snr_m=[0 snr_m 1e100];
[Mat_R2,Mat_snr]=meshgrid(R_set,snr_m');
sigma=sqrt((1-rho).*snr./2);
eR_m = normrnd(0,sigma,K,N);
eI_m = normrnd(0,sigma,K,N);
heR_m = normrnd(0,sqrt(snr./2),K,N);
heI_m = normrnd(0,sqrt(snr./2),K,N);

for ind_lambda=1:length(lambda_set)
    lambda=lambda_set(ind_lambda);
    
    pdf_snr=@(g) exp(-g./snr)./snr;
    pdf_cond=@(g,ge) exp(-(g+rho.*ge-2*sqrt(rho.*g.*ge))./((1-rho).*snr)).*besseli(0,2.*sqrt(rho.*g.*ge)./((1-rho).*snr),1)./((1-rho).*snr);
    
    V_2=@(R,n,ge2) R.*(R>=n).*marcumq(sqrt(rho*2.*ge2./((1-rho).*snr)),sqrt(2.*(2.^R-1)./((1-rho).*snr))).*(n>0)-lambda.*(n>0).*(R>=n);
    V_2c=@(R,n,ge2,g1) (R+log2(1+g1)).*(R>=n).*marcumq(sqrt(rho*2.*ge2./((1-rho).*snr)),sqrt(2.*(2.^R-1)./((1-rho).*snr))).*(n>0)-lambda.*(n>0).*(R>=n);
    MParity_Rate2=zeros(length(R_set),length(x)+2);
    tic
    for i=1:length(R_set)
        %  tic
        [val0,ind]=max((V_2(Mat_R,R_set(i),Mat_g))');
        Parity_Rate2(i,:)=R_set(ind);
        for j=1:length(x)
            val=V_2c(R_set(ind)',R_set(i),x,x(j));
            xtemp=val.*pdf_snr(x);
            V_av2(i,j)=sum(xtemp.*w)*Dx/2;
        end
        %  toc
    end
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
            val_temp=val_temp+V_1(R_set(j),x(i));
            temp(j)=sum(val_temp.*w)*th_d/2;
        end
        [val,ind]=max(temp);
        V{1}(i)=val;
        Parity_Rate1(i)=R_set(ind);
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




% xtemp=squeeze(V{1})'.*pdf_snr(x);
% V_av1=sum(xtemp.*w)*Dx/2;
% toc
% test_lambda=-lambda+V_av1

1













