function p=proba_joint_mutual_inform_1(SNR,a,b,c,d)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%Cette fonction calcul la probabilite 
%Pr[a<log2(1+X)<b and c<log2(1+X)<d]
%avec un canal rayleigh fading E[X]=snr
%SNR en dB
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

snr=10^(SNR/10);

a1=max(a,c);
b1=min(b,d);

if(a1>b1)
    p=0;
    return
else
    p=exp(-max((2^a1-1),0)/snr)-exp(-max((2^b1-1),0)/snr);
    return
end