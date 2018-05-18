function out=WEP_parfait(snr,rate,g)

% out=snr<=log2(rate);
load MI_16QAM.mat;
MI_16QAM=[0 MI_16QAM(1:end-1)' 4];
th=[0 10.^((SNR(1:end-1))'./10) 1e10];
Qam_th=@(R) interp1(MI_16QAM,th,R);

snr_th=Qam_th(rate);
out=exp(-g.*(snr./snr_th-1));
out(find(snr-snr_th<0))=1;   
out(find(snr_th==0))=0; 


end