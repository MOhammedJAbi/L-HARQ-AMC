 clear;
clc;clf;

load WEP_Turbo_16QAM_rate1.5poincenage6_3gpp.mat
%[val,nan_v]=max(isnan(WEP_set(:,1)));
%[val,nan_v]=min((WEP_set(:,1))>=1e-4);
%j=nan_v;
%j=(nan_v-1)*(nan_v>1)+length(SNRLin)*(nan_v==1);
% SNRLin_m{i}=[0 SNRLin(1:j) 1e100];
% WEP_m{i}=[1 (WEP_set((1:j),1))' 0];
i=1;
[C,IA,IC] = unique(WEP_set(:,1), 'last') ;
output_wep = WEP_set((sort(IA)),1);
output_SNRLin=SNRLin(sort(IA));
[val,nan_v]=max(isnan(output_wep));
j=(nan_v-1)*(nan_v>1)+length(IA)*(nan_v==1);
SNRLin_m{i}=[1e-50 output_SNRLin(1:j) 1e100];
WEP_m{i}=[1 output_wep(1:j)' 1e-50];

load WEP_Turbo_16QAM_rate2.25poincenage6_3gpp.mat
i=2;
[C,IA,IC] = unique(WEP_set(:,1), 'last') ;
output_wep = WEP_set((sort(IA)),1);
output_SNRLin=SNRLin(sort(IA));
[val,nan_v]=max(isnan(output_wep));
j=(nan_v-1)*(nan_v>1)+length(IA)*(nan_v==1);
SNRLin_m{i}=[1e-50 output_SNRLin(2:j-1) 1e100];
WEP_m{i}=[1 output_wep(2:j-1)' 1e-50];


load WEP_Turbo_16QAM_rate3poincenage6_3gpp.mat
i=3;
[val,nan_v]=min((WEP_set(:,1))>=1e-4);
[C,IA,IC] = unique(WEP_set(:,1), 'last') ;
output_wep = WEP_set((sort(IA)),1);
output_SNRLin=SNRLin(sort(IA));
[val,nan_v]=max(isnan(output_wep));
j=(nan_v-1)*(nan_v>1)+length(IA)*(nan_v==1);
SNRLin_m{i}=[1e-50 output_SNRLin(1:j-1) 1e100];
WEP_m{i}=[1 output_wep(1:j-1)' 1e-50];
% semilogy(10.*log10(SNRLin_m{i}),WEP_m{i},'r',10.*log10(SNRLin),WEP_set(:,1),'-.b')
%  axis([-4 30 1e-5 1])

load WEP_Turbo_16QAM_rate3.75poincenage6_3gpp.mat
i=4;
[C,IA,IC] = unique(WEP_set(:,1), 'last') ;
output_wep = WEP_set((sort(IA)),1);
output_SNRLin=SNRLin(sort(IA));
[val,nan_v]=max(isnan(output_wep));
j=(nan_v-1)*(nan_v>1)+length(IA)*(nan_v==1);
SNRLin_m{i}=[1e-50 output_SNRLin(2:j-1) 1e100];
WEP_m{i}=[1 output_wep(2:j-1)' 1e-50];


load AMC_policy_16QAM_Turbo_code_dop_0.05rho_0.95156.mat

SNRdB=[-5:30];
% correlation coefficient
dop=0.05
rho=(besselj(0,2*pi*dop)).^2

Rate_set=[1.5 2.25 3 3.75];

K=1;
N=1e5;
[x0, w0] = GaussLegendre(200);

%% Computation of throughput
for ind_snr=1:length(SNRdB)
    ind_snr
    snr=10.^(SNRdB(ind_snr)./10);
    sigma=sqrt((1-rho).*snr./2);
    Dx=-log(10^-15)*snr;
    gammae=(x0+1)*Dx/2;
 
    phi = unifrnd(0,2.*pi,K,N);
    eR = normrnd(0,sigma,K,N);
    eI = normrnd(0,sigma,K,N);
    
    
%     for ind_SNRhat=1:length(gammae)
%         gamma=rho.*gammae(ind_SNRhat)+eI.^2+eR.^2-2.*sqrt(rho).*eR.*sqrt(gammae(ind_SNRhat)).*cos(phi)-2.*sqrt(rho).*eI.*sqrt(gammae(ind_SNRhat)).*sin(phi);
%         for j=1:length(Rate_set)
%             wep_v=interp1(10.*log10(squeeze(SNRLin_m{j})),10.*log10(squeeze(WEP_m{j})),10.*log10(gamma));
%             wep_v=10.^(wep_v./10);
%             th(j)=Rate_set(j).*sum(binornd(1, 1-wep_v))./N;
%         end
%         [val,ii]=max(th);
%         policy_tem(ind_SNRhat)=Rate_set(ii);
%     end
%    
%     
%     policy(ind_snr,:)=[policy_tem(1) policy_tem Rate_set(end)];
%     snr_range(ind_snr,:)=[1e-50 gammae' 1e100];
    
    tt=1;
    eR = normrnd(0,sigma,K,N);
    eI = normrnd(0,sigma,K,N);
    heR = normrnd(0,sqrt(snr./2),K,N);
    heI =normrnd(0,sqrt(snr./2),K,N);
    gammae=heR.^2+ heI.^2;
    gamma=(sqrt(rho).*heR-eR).^2+(sqrt(rho).*heI-eI).^2;
    Optimal_rate=interp1(snr_range(ind_snr,:),policy(ind_snr,:),gammae,'nearest');
    wep=zeros(1,N);
    I1=find(Optimal_rate==1.5);
    wep(I1)=interp1(10.*log10(squeeze(SNRLin_m{1})),10.*log10(squeeze(WEP_m{1})),10.*log10(gamma(I1)));
    I2=find(Optimal_rate==2.25);
    wep(I2)=interp1(10.*log10(squeeze(SNRLin_m{2})),10.*log10(squeeze(WEP_m{2})),10.*log10(gamma(I2)));
    I3=find(Optimal_rate==3);
    wep(I3)=interp1(10.*log10(squeeze(SNRLin_m{3})),10.*log10(squeeze(WEP_m{3})),10.*log10(gamma(I3)));
    I4=find(Optimal_rate==3.75);
    wep(I4)=interp1(10.*log10(squeeze(SNRLin_m{4})),10.*log10(squeeze(WEP_m{4})),10.*log10(gamma(I4)));
    wep=10.^(wep./10);
    throughput_amc(ind_snr)=sum(Optimal_rate.*binornd(1, 1-wep))./N;
    

end

%save(['AMC_policy_16QAM_Turbo_code_dop_',num2str(dop),'rho_',num2str(rho),'2.mat'],'throughput_amc', 'snr_range','policy','SNRdB','dop','r')

