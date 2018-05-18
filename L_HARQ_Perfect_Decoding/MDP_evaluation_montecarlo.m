function [T,ac_ret,ac_sc,ac_drop,ac_no_need]=MDP_evaluation_montecarlo(R,Policy,SNR,M)


N=100000; % nombre de realisation de canal
snr=10^(SNR/10);
X=exprnd(snr,1,N);
%step=length(Policy);
%step=step-1;
ac_ret=0;
ac_sc=0;
ac_no_need=0;
ac_drop=0;
Rew=0;
I_s=zeros(1,M+1);

for i=1:M+1
    I_s(i)=(i-1)*R/M;
end

i=2;

I=log2(1+X(1));

while(i<N)
    
    K=find(I>=I_s);
    %I_s(K(length(K)))
    ind=K(end);
    alpha=Policy(ind);
    if(alpha~=2)
        ac_no_need=ac_no_need+(I>=R);
        ac_drop=ac_drop+(I<R&&alpha==0);
        ac_sc=ac_sc+(I<R&&alpha~=0);
        if(I+alpha*log2(1+X(i))>=R)
            Rew=Rew+R;
            I=(1-alpha)*log2(1+X(i));
            i=i+1;
        else
            I=log2(1+(1-alpha)*X(i)/(1+alpha*X(i)));
            i=i+1;
        end
    else
        ac_ret=ac_ret+1;
        if(I+log2(1+X(i))>=R)
            Rew=Rew+R;
            I=log2(1+X(i+1));
            i=i+2;
        else
            I=log2(1+X(i+1));
            i=i+2;
        end
    end
        
end

T=Rew/N;
KK=ac_ret+ac_sc+ac_drop+ac_no_need;
ac_ret=ac_ret/KK;
ac_sc=ac_sc/KK;
ac_no_need=ac_no_need/KK;
ac_drop=ac_drop/KK;

% 
% lmont = 1e7;
% 
% 
% gbar = 10^(gdB/10);
% 
% 
% qq = rho(1,:)==0;
% if sum(qq) == K
%     Po = 1;
%     Eta = 0;
%     return
% end
% 
% for i = 1:10
%     
% cSD = random('exp',gbar,lmont,K);
% Id = rho(1,1)*log2(1+cSD(:,1));
% totR = ones(lmont,1).*rho(1,1);
% for k = 2:K
%     cId = Id<1;
%     r = interp1(x,rho(:,k),Id.*cId).*cId;
%     Id = Id + r.*log2(1+cSD(:,k));
%     totR = totR + r;        
% end
% den(i) = mean(totR);
% num(i) = mean(Id>=1);
% po(i) = mean((Id<1));
% 
% end
% Po = mean(po);
% Eta = mean(num)/mean(den);