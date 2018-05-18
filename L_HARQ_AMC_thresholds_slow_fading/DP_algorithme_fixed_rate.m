function [Throughput,Parity_Rate,gamma]=DP_algorithme_fxed_rate(R_set,snrth_set,gamma_op,snr,a,K,L)


%% Dynamic programming approach

% SNR pdf, av_snr=channel average SNR
pdf_snr=@(x,av_snr) exp(-x./av_snr)./av_snr;
% SNR cdf
cdf_snr=@(x,av_snr) 1-exp(-x./av_snr);
% PER average
PER=@(l,x,av_snr) pdf_snr(x,av_snr).*WEP(x,snrth_set(l),a);
% Legend integral 
lint =200;
[x0, w] = GaussLegendre(lint);
% decoding threshold, used in PER expression 
snrth=@(R) 2.^R-1;
% gamma
gamma=[0 logspace(log10(min(R_set(1),snrth_set(1))),log10(gamma_op(end)),1e3)];
% throughput computation
l=1
Throughput(l,:)=R_set(l).*(1-WEP(gamma,snrth_set(l),a));
V=@(l,gamma,Rp) (R_set(l)-Rp).*(WEP(gamma,snrth_set(l),a)-WEP(gamma,snrth(R_set(l)-Rp),a));
for l=2:L
    Rp_set=0;
    Rp_set=[0 R_set(1:l)];
    I=gamma';
    [Mat_R,Mat_I]=meshgrid(Rp_set,I);
    [val,ind]=max((V(l,Mat_I,Mat_R))');
   
    Valu=val'.*(1-WEP(I,snrth_set(l),a))+R_set(l).*(1-WEP(I,snrth_set(l),a))...
        +R_set(l).*(1-WEP(I,snrth_set(l),a)).*WEP(I,snrth_set(l),a);
    Throughput(l,:)=(Valu./(1+WEP(I,snrth_set(l),a)));
    Parity_Rate(l,:)=Rp_set(ind);
end
    

    
    
