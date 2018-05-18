function [lambda_k,xop,counter,Pout]=Policy_iteration_algorithme(R_set,snr,M,Rp_set,alpha0,snr_max,rho,ind_hp,Max_iter,s) %(R_set,snrth_set,snr,M,Rp_set,Rp0,snr_max,rho,ind_hp,Max_iter);

%function [lambda_k,xop,counter]=Policy_iteration_algorithme(R1,R2,SNR1,SNR2,step,beta1,beta2)

%%%%%%%%% Policy iteration %%%%%%%%

 
LR=length(R_set);
Lp=length(Rp_set);
alpha0((M+1))=R_set(L);


%Prob_trans_alpha_set=zeros(s,s,L);
%tic
P_k=Matrix_transition(alpha0,M,snr,snr_max,s,rho,ind_hp);
%tt=toc
% aaa=sum(P_k,2);
% bb=[1:s];
% bbb=[aaa alpha0' bb']
%  max(bbb(:,1))
%P_k(10,:)
% aaa(5:20)
rew_k=zeros(s,1);

for i=1:s
    rew_k(i)=reward_expectation(i,alpha0(i),R,SNR,M);
end 

I=eye(size(P_k));
K=inv(I-P_k);
c=K*rew_k;% c=lambda^(-1)+h
lambda_k=c(ind_hp)/(sum(K(ind_hp,:)));
h_k=c-K*lambda_k*ones(s,1);


counter=0;
opt_matrix=ones(s,L);
Prob_trans_alpha_set=zeros(s,s,L);
Reward_alpha_set=zeros(s,L);
for j=1:L
    for i=1:s
        Reward_alpha_set(i,j)=reward_expectation(i,alpha_set(j),R,SNR,M);
    end
    
    %tic
    Prob_trans_alpha_set(:,:,j)=Matrix_transition(alpha_set(j)*ones(s),M,R,SNR,ind_hp);% probability_trans_sc_adp(i,alpha_set(j),R,SNR,step);
    %ttt=toc
end
    
while(counter<Max_iter)
    for i=1:s
        for j=1:L
            a=Reward_alpha_set(i,j);
            b=Prob_trans_alpha_set(i,:,j);
            opt_matrix(i,j)=a+b*h_k;
        end
        l=find(opt_matrix(i,:)==max(opt_matrix(i,:)));
        xop(i)=alpha_set(l(end));
    end
    xop((M+1):(M+1):(M+1)*(M+2))=0;
    xop((M+1)*(M+2):end)=0;
    %%%%%%%%%%%%%%%%  h_k_1    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %xop=ones(1,N12);
    
    P_k_1=Matrix_transition(xop,M,R,SNR,ind_hp); %Matrix_transition(A1,xop,R1,R2,SNR1,SNR2,step);
   
    
    rew_k_1=zeros(s,1);
    for i=1:s
        rew_k_1(i)=reward_expectation(i,xop(i),R,SNR,M); %(i,A1(i),xop(i),R1,R2,SNR1,SNR2,step,beta1,beta2);
    end
    I=eye(size(P_k_1));
%    A=I-P_k_1;
%     si=size(A)
%     ran=rank(A)
%     de=det(A) 
    K=inv(I-P_k_1);
    c=K*rew_k_1;  % c=lambda^(-1)+h
    lambda_k_1=c(ind_hp)/(sum(K(ind_hp,:)));
    h_k_1=c-lambda_k_1*K*ones(s,1);
    %aa=h_k-h_k_1;
    %max(abs(aa))
    if(isequal(h_k,h_k_1)&&lambda_k==lambda_k_1)
        break
    end
    h_k=h_k_1;
    lambda_k=lambda_k_1;
    counter=counter+1;
end


%%%%%%%%%%%   outage   %%%%%%%%%%%%%%%%%%%

rew_of_re_k=zeros(s,1);

for i=1:s
    rew_of_re_k(i)=Np_expectation(i,xop(i),R,M,SNR);
end 

O=K*rew_of_re_k;
zeta_p=O(ind_hp)/(sum(K(ind_hp,:)));
Pout=1-lambda_k/(R*zeta_p);
 
