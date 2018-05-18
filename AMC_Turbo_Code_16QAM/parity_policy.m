function [parity_policy,gamma_interval] = parity_policy(epsilon_Set)

% Legend integral
number_of_point= 4e2;
[x0, w] = GaussLegendre(number_of_point);


%clf;
SNRLin=0;
WEP_set=0;

%% R=1.5
load WEP_Turbo_16QAM_rate1.5poincenage6_3gpp.mat
SNRLin(SNRLin==0)=1e-50;
WEP_set(isnan(WEP_set))=1e-50;%min(WEP_set(WEP_set>0));
WEP_set(WEP_set==0)=1e-50;
[Rate_parity_Set_M,snr_M]=meshgrid(Rate_parity_Set,SNRLin);
FInterpolation1 = scatteredInterpolant(10.*log10(snr_M(:)),Rate_parity_Set_M(:),10.*log10(WEP_set(:)));
WEP1=@(gamma,rate_parity) 10.^(FInterpolation1(10.*log10(gamma),rate_parity)./10);
%% optimization of epsilon
Dxmin=SNRLin(2);
Dxmax=SNRLin(end-1);
IM1=[SNRLin(1) ((x0*(Dxmax-Dxmin)+(Dxmax+Dxmin))./2)' SNRLin(end)];
[Rate_parity_Set_M,I_M]=meshgrid(Rate_parity_Set,IM1);
WEP_set_M = WEP1(I_M,Rate_parity_Set_M) ; %  F(I_M,Rate_parity_Set_M);
M=length(IM1);
L=length(Rate_parity_Set);
WEP_set_M_N = WEP_set_M./WEP1(I_M,zeros(M,L));


for t=1:length(epsilon_Set)
    
    epsilon=epsilon_Set(t);
    V_temp_M=WEP_set_M_N-epsilon<=0;
    [m,ind]=max(V_temp_M>0,[],2);
    parity_epsilon1(t,:)=Rate_parity_Set(ind);
    parity_policy(1,t,:)=Rate_parity_Set(ind);
end
gamma_interval(1,:)=IM1;
%parity_epsilon1(parity_epsilon1==0)=min(Rate_parity_Set(Rate_parity_Set>0));% Rate_parity_Set
% t=30;
% epsilon_Set(t)
% loglog(IM1(2:end-1),WEP1(IM1(2:end-1),parity_epsilon1(t,(2:end-1)))./ WEP1(IM1(2:end-1),zeros(1,M-2)) )


%% R=2.25
SNRLin=0;
WEP_set=0;
load WEP_Turbo_16QAM_rate2.25poincenage6_3gpp.mat
SNRLin(SNRLin==0)=1e-50;
WEP_set(isnan(WEP_set))=1e-50;
WEP_set(WEP_set==0)=1e-50;
[Rate_parity_Set_M,snr_M]=meshgrid(Rate_parity_Set,SNRLin);
FInterpolation2 = scatteredInterpolant(10.*log10(snr_M(:)),Rate_parity_Set_M(:),10.*log10(WEP_set(:)));
WEP2=@(gamma,rate_parity) 10.^(FInterpolation2(10.*log10(gamma),rate_parity)./10);
Dxmin=SNRLin(2);
Dxmax=SNRLin(end-1);
IM2=[SNRLin(1) ((x0*(Dxmax-Dxmin)+(Dxmax+Dxmin))./2)' SNRLin(end)];
[Rate_parity_Set_M,I_M]=meshgrid(Rate_parity_Set,IM2);
WEP_set_M = WEP2(I_M,Rate_parity_Set_M) ; %  F(I_M,Rate_parity_Set_M);
M=length(IM2);
L=length(Rate_parity_Set);
WEP_set_M_N = WEP_set_M./WEP2(I_M,zeros(M,L));

for t=1:length(epsilon_Set)
    
    epsilon=epsilon_Set(t);
    V_temp_M=WEP_set_M_N-epsilon<=0;
    [m,ind]=max(V_temp_M>0,[],2);
    parity_epsilon2(t,:)=Rate_parity_Set(ind);
    parity_policy(2,t,:)=Rate_parity_Set(ind);
end
gamma_interval(2,:)=IM2;

%% R=3
SNRLin=0;
WEP_set=0;
load WEP_Turbo_16QAM_rate3poincenage6_3gpp.mat
SNRLin(SNRLin==0)=1e-50;
WEP_set(isnan(WEP_set))=1e-50;
WEP_set(WEP_set==0)=1e-50;
[Rate_parity_Set_M,snr_M]=meshgrid(Rate_parity_Set,SNRLin);
FInterpolation3 = scatteredInterpolant(10.*log10(snr_M(:)),Rate_parity_Set_M(:),10.*log10(WEP_set(:)));
WEP3=@(gamma,rate_parity) 10.^(FInterpolation3(10.*log10(gamma),rate_parity)./10);
Dxmin=SNRLin(2);
Dxmax=SNRLin(end-1);
IM3=[SNRLin(1) ((x0*(Dxmax-Dxmin)+(Dxmax+Dxmin))./2)' SNRLin(end)];
[Rate_parity_Set_M,I_M]=meshgrid(Rate_parity_Set,IM3);
WEP_set_M = WEP3(I_M,Rate_parity_Set_M) ; %  F(I_M,Rate_parity_Set_M);
M=length(IM3);
L=length(Rate_parity_Set);
WEP_set_M_N = WEP_set_M./WEP3(I_M,zeros(M,L));

for t=1:length(epsilon_Set)
    epsilon=epsilon_Set(t);
    V_temp_M=WEP_set_M_N-epsilon<=0;
    [m,ind]=max(V_temp_M>0,[],2);
    parity_epsilon3(t,:)=Rate_parity_Set(ind);
    parity_policy(3,t,:)=Rate_parity_Set(ind);
end
gamma_interval(3,:)=IM3;



%% R=3.75
SNRLin=0;
WEP_set=0;
load WEP_Turbo_16QAM_rate3.75poincenage6_3gpp.mat
SNRLin(SNRLin==0)=1e-50;
WEP_set(isnan(WEP_set))=1e-50;
WEP_set(WEP_set==0)=1e-50;
[Rate_parity_Set_M,snr_M]=meshgrid(Rate_parity_Set,SNRLin);
FInterpolation4 = scatteredInterpolant(10.*log10(snr_M(:)),Rate_parity_Set_M(:),10.*log10(WEP_set(:)));
WEP4=@(gamma,rate_parity) 10.^(FInterpolation4(10.*log10(gamma),rate_parity)./10);
Dxmin=SNRLin(2);
Dxmax=SNRLin(end-1);
IM4=[SNRLin(1) ((x0*(Dxmax-Dxmin)+(Dxmax+Dxmin))./2)' SNRLin(end)];
[Rate_parity_Set_M,I_M]=meshgrid(Rate_parity_Set,IM4);
WEP_set_M = WEP4(I_M,Rate_parity_Set_M) ; %  F(I_M,Rate_parity_Set_M);
M=length(IM4);
L=length(Rate_parity_Set);
WEP_set_M_N = WEP_set_M./WEP4(I_M,zeros(M,L));

for t=1:length(epsilon_Set)
    
    epsilon=epsilon_Set(t);
    V_temp_M=WEP_set_M_N-epsilon<=0;
    [m,ind]=max(V_temp_M>0,[],2);
    parity_epsilon4(t,:)=Rate_parity_Set(ind);
    parity_policy(4,t,:)=Rate_parity_Set(ind);
end
gamma_interval(4,:)=IM4;

