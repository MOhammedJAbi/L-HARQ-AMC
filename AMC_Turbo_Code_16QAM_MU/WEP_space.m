function [FInterpolation1,FInterpolation2,FInterpolation3,FInterpolation4] = WEP_space(gamma)

gamma=gamma+1;
% This function estimates the WEP

SNRLin=0;
WEP_set=0;
load WEP_Turbo_16QAM_rate1.5poincenage6_3gpp.mat
SNRLin(SNRLin==0)=1e-50;
WEP_set(isnan(WEP_set))=1e-50;
WEP_set(WEP_set==0)=1e-50;
[Rate_parity_Set_M,snr_M]=meshgrid(Rate_parity_Set,SNRLin);
FInterpolation1 = scatteredInterpolant(10.*log10(snr_M(:)),Rate_parity_Set_M(:),10.*log10(WEP_set(:)));

SNRLin=0;
WEP_set=0;
load WEP_Turbo_16QAM_rate2.25poincenage6_3gpp.mat
SNRLin(SNRLin==0)=1e-50;
WEP_set(isnan(WEP_set))=1e-50;
WEP_set(WEP_set==0)=1e-50;
[Rate_parity_Set_M,snr_M]=meshgrid(Rate_parity_Set,SNRLin);
FInterpolation2 = scatteredInterpolant(10.*log10(snr_M(:)),Rate_parity_Set_M(:),10.*log10(WEP_set(:)));

SNRLin=0;
WEP_set=0;
load WEP_Turbo_16QAM_rate3poincenage6_3gpp.mat
SNRLin(SNRLin==0)=1e-50;
WEP_set(isnan(WEP_set))=1e-50;
WEP_set(WEP_set==0)=1e-50;
[Rate_parity_Set_M,snr_M]=meshgrid(Rate_parity_Set,SNRLin);
FInterpolation3 = scatteredInterpolant(10.*log10(snr_M(:)),Rate_parity_Set_M(:),10.*log10(WEP_set(:)));


SNRLin=0;
WEP_set=0;
load WEP_Turbo_16QAM_rate3.75poincenage6_3gpp.mat
SNRLin(SNRLin==0)=1e-50;
WEP_set(isnan(WEP_set))=1e-50;
WEP_set(WEP_set==0)=1e-50;
[Rate_parity_Set_M,snr_M]=meshgrid(Rate_parity_Set,SNRLin);
FInterpolation4 = scatteredInterpolant(10.*log10(snr_M(:)),Rate_parity_Set_M(:),10.*log10(WEP_set(:)));








