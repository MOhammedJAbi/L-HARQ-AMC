clear;
clc;clf;


% correlation coefficient

pdf_snr=@(x,av_snr) exp(-x./av_snr)./av_snr;
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
SNR=[-5:40];


%% Computation of throughput
for ind_snr=1:length(SNR)
    ind_snr
    snr=10.^(SNR(ind_snr)./10);
    
    Dx=snr.*8;
    x = (x0+1)*Dx/2;
    b=pdf_snr(x,snr).*Qam_r(x);
    Er_MI_64QAM_fading(ind_snr)=sum(b.*w)*Dx/2;

    
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

save Er_MI_64QAM.mat SNR Er_MI_64QAM_fading

plot(SNR,Er_MI_64QAM_fading,'r')
hold on
load Er_MI_16QAM.mat

plot(SNR,Er_MI_16QAM_fading,'-.b')