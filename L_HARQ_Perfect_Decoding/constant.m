clear;
clc;
K=8;
exp_cont_outage=6;
cont_outage=10^-exp_cont_outage;
rate=4 %[0.5 1 1.5 3 4.5];
m=1;
N=1e4;
cdf=@(P,snr,x) 1-gammainc(m.*(2.^x-1)./(snr.*P),m,'upper');
pdf=@(P,x,snr) m^m.*log(2).*2.^x.*(2.^x-1).^(m-1).*exp(-m.*(2.^x-1)./(snr.*P))./((snr.*P).^m.*gamma(m));
%snrdB=[-5:0.1:40]; %version 1
snrdB=[-5:1:40]; %version 2
for t=1:length(rate)
    R=rate(t);
    pas=R/N;
    x=linspace(0,R-pas,N)+pas/2;
    for i=1:length(snrdB)
        snr=10.^(snrdB(i)/10);
        time=1;
        for tt=1:K
            if tt==1
                f(tt,i,t)=cdf(1,snr,R);
            else
                tem_fft=1;
                for l=1:tt
                    tem_fft=fft(pdf(1,x,snr),2^ceil(log2(tt*length(x)-1))).*tem_fft;
                end
                FF=ifft(tem_fft);
                f(tt,i,t)=(pas^tt)*sum(FF(1:length(x)));
            end
            throughput(tt,i,t)=rate(t)*(1-f(tt,i,t))./time;
            if f(tt,i,t)>cont_outage
                throughput_cont(tt,i,t)=0;
            else
                throughput_cont(tt,i,t)=throughput(tt,i,t);
            end
            time=time+f(tt,i,t);
        end
    end
    
end
 
save('constant_rate')