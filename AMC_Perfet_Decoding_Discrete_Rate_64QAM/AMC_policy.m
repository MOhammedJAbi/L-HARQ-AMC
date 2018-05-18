 clear;
clc;clf;


% correlation coefficient
dop=0.1
rho=(besselj(0,2*pi*dop)).^2

% pdf joint, av_snr=channel average SNR
pdf_joint=@(x,y,av_snr) exp(-(x+y-2*sqrt(rho.*x.*y))./((1-rho).*av_snr)).*besseli(0,2.*sqrt(rho.*x.*y)./((1-rho).*av_snr),1)./((1-rho).*av_snr.^2);
pdf_cond=@(x,y,av_snr) (1-marcumq(sqrt(rho*2.*x./((1-rho).*av_snr)),sqrt(2.*y./((1-rho).*av_snr)))).*exp(-x./av_snr)./av_snr;
pdf_snr=@(x,av_snr) exp(-x./av_snr)./av_snr;
% SNR cdf
cdf_snr=@(x,av_snr) 1-exp(-x./av_snr);
% optimal threshold
load MI_64QAM_SNRdB=-20_30_sigmap2=0_0.mat
MI_16QAM=[0 MI_PN_BLT(1:end-5)' 6];
th=[0 10.^((SNRdB(1:end-5))./10) 1e10];
Qam_th=@(R) interp1(MI_16QAM,th,R)
Qam_r=@(gamma) interp1(th,MI_16QAM,gamma)

 
f_threshold1=@(x,R,av_snr) R.*marcumq(sqrt(rho*2.*x./((1-rho).*av_snr)),sqrt(2.*Qam_th(R)./((1-rho).*av_snr)));

% Legend integral 
lint = 2e2;
[x0, w] = GaussLegendre(lint);
SNRdB=[-5:40];


%% Computation of throughput
for ind_snr=1:length(SNRdB)
    ind_snr
    snr=10.^(SNRdB(ind_snr)./10);
    
    max_rate=4;
    Dx=snr.*8;
    x = (x0+1)*Dx/2;
    y=[2.25;3;3.75;4.5;5.25;6];%(linspace(0.01,max_rate-1e-5,1e2))'
    

     [xmat,ymat] = meshgrid(x,y);
     [f,I]=max(f_threshold1(xmat,ymat,snr));
    policy(ind_snr,:)=y(I);
    snr_range(ind_snr,:)=x;
    b=pdf_snr(x,snr);
    throughput(ind_snr)=sum(f'.*b.*w)*Dx/2;

    
    %toc
%     
%     [f,I]=max(f_threshold1(xx,y,snr))
%          figure(2)
%     plot(y,f_threshold1(xx,y,snr),'r',y,f_threshold2(xx,y,snr),'-.b')
% tic
%     [xmat,ymat] = meshgrid(x,y);
%     
%     f_threshold1(xmat,ymat,snr);
%     [f,I]=max(f_threshold1(xmat,ymat,snr));
%     toc
%     xx=x(1)
%     figure(1)
%     plot(y,f_threshold(xx,y,snr))
%     figure(2)
%     plot(y,f_threshold1(xx,y,snr))
    
%finding optimal threshold   

% 

    
end

save(['AMC_policy_64QAM_Perfect_Decoding_Discrete_R_dop_',num2str(dop),'rho_',num2str(rho),'.mat'])
% figure(1)

% gamma_op=10*log10(gamma_op);
% plot(SNRdB,gamma_op(:,2),'r')%SNRdB,gamma_op_app(:,l),'-.b')
% hold on
% plot(SNRdB,gamma_op(:,3),'g')
% plot(SNRdB,gamma_op(:,4),'b')
% plot(SNRdB,gamma_op(:,5),'c')

grid
figure(2)

plot(SNRdB,throughput,'-r')

hold on

load Er_MI_16QAM.mat

plot(SNR,Er_MI_16QAM_fading,'-.b')

grid

figure(1)
i=1
plot(snr_range(i,:),policy(i,:),'r')
hold on 
plot(snr_range(i,:),log2(1+snr_range(i,:)),'b')
hold on
i=10
plot(snr_range(i,:),policy(i,:),'b')
hold on
i=20
plot(snr_range(i,:),policy(i,:),'c')
hold on
i=30
plot(snr_range(i,:),policy(i,:),'m')
hold on
i=36
plot(snr_range(i,:),policy(i,:),'g')







% 
% figure(2)
% plot(SNRdB,T_corre,'r')
% grid
% 
% figure(2)
% semilogy(SNRdB, outage(:,1),SNRdB, outage(:,2),SNRdB, outage(:,3),SNRdB, outage(:,4),SNRdB, outage(:,5))
% legend('l=1','l=2','l=3','l=4','l=5')







