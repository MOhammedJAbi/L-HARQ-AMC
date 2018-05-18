function P=probability_trans_sc_adp(i,alpha_i,M,R,SNR)


snr=10^(SNR/10);
% %%%%%%%%%%%%%%%%%%%%%%%       states       %%%%%%%%%%%%%%%%%%%%%%
% 
% 
% I=zeros(1,step+1);
% 
% for i=1:step+1
%     I(i)=(i-1)*R/step;
% end

%%%%%%%%%%%%%%%%%%%%%%%     transition matrix     %%%%%%%%%%%%%%%%%%%

P=zeros(1,M);
for i=1:M
    P(i)=
    
% if (i~=2*(step+1)+1)
% 
% P1=zeros(1,step+1);
% P2=zeros(1,step+1);
% 
% 
% 
% 
%     for j=1:step+1
%         P1(j)=prb_trans_sc_adpt(i,j-1,1,alpha_i,R,SNR,step);
%         
%         P2(j)=prb_trans_sc_adpt(i,j-1,0,alpha_i,R,SNR,step);
%     end
% p1=prb_trans_sc_adpt(i,2*(step+1)+1,1,alpha_i,R,SNR,step);
% p2=prb_trans_sc_adpt(i,2*(step+1)+1,0,alpha_i,R,SNR,step);
% P=[P1 P2 p1 p2];
% return
% end
% 
% 
% if (i==2*(step+1)+1)
%     P=zeros(1,2*(step+1)+2);
%     
%     for i=1:step+1
%         P(i+step+1)=prb_trans_sc_adpt(2*step+3,i-1,0,1,R,SNR,step);
%     end
% 
%     P(2*(step+1)+1)=exp(-(2^R-1)/snr);
% end
return

