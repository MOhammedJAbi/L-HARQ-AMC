clc
clear

% rate set
pas=0.75
R_set=[pas:pas:3.75];
L=length(R_set);
% decoding threshold, used in PER expression 
snrth=@(R) 2.^R-1;
snrth_set=snrth(R_set)
% systematic rates
Rp_set=[0 R_set];

SNRdB=[-5:30];
K=2;  %nombre maximal de transmission
M=31;  %discritisation de l'information mutuelle
s=(M+1)*(L*M+1);  %nombre total d'etat
ind_hp=1;  %indice de hp
Max_iter=15; %nombre maximal d'iteration

% correlation coefficient
rho=0.8 %(besselj(0,2*pi*1))
% definition of snr_max
epsilon=1e-6;
f_threshold=@(x,k,av_snr) (1-marcumq(sqrt(rho*2.*x./((1-rho).*av_snr)),sqrt(2.*snrth_set(k)./((1-rho).*av_snr))))./epsilon-1;
gamma_max=snrth_set(L)*3e2;

T_corre=zeros(1,length(SNRdB));
conter_op=zeros(1,length(SNRdB));
Policy_corre=zeros(length(SNRdB),s);
Rp0=zeros(1,s);
%alpha0(10)=alpha_step(1);

%alpha0=zeros(1,s);
%alpha0=alpha_set(randi(length(alpha_set),s,1));

for j=1:length(SNRdB)
    j
    snr=10.^(SNRdB(j)./10);
    snr_max(j)=bisection(@(x) f_threshold(x,L,snr), 0,gamma_max,100, 1e-3, 1e-2);
    
    [lambda_k,xop,counter]=Policy_iteration_algorithme(R_set,snrth_set,snr,M,Rp_set,Rp0,snr_max,rho,ind_hp,Max_iter,s);
    T_corre(j)=lambda_k;
    Policy_corre(j,:)=xop;
    conter_op(j)=counter;
    Rp0=xop;
end

plot(SNRdB,snr_max)


%save('fixed_rate_var_snr_R_V8')
   
%save(['TS_fixed_rate_var_snr_R_', num2str(R),'_K_', num2str(K),'_MI_', num2str(M),'_Alpha_', num2str(L),'.mat'])