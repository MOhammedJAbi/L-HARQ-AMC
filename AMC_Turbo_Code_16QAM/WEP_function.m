function [ WEP] = WEP_function(transmission_rate,parity_rate,gamma)

% This function estimates the WEP

if transmission_rate==1.5
    
    load WEP_Turbo_16QAM_rate1.5poincenage6_3gpp.mat
    SNRLin(SNRLin==0)=1e-50;
    WEP_set(isnan(WEP_set))=1e-50;
    WEP_set(WEP_set==0)=1e-50;
    [Rate_parity_Set_M,snr_M]=meshgrid(Rate_parity_Set,SNRLin);
    FInterpolation1 = scatteredInterpolant(10.*log10(snr_M(:)),Rate_parity_Set_M(:),10.*log10(WEP_set(:)));
    WEP= 10.^(FInterpolation1(10.*log10(gamma),parity_rate)./10);
    return

elseif transmission_rate==2.25
    
    load WEP_Turbo_16QAM_rate2.25poincenage6_3gpp.mat
    SNRLin(SNRLin==0)=1e-50;
    WEP_set(isnan(WEP_set))=1e-50;
    WEP_set(WEP_set==0)=1e-50;
    [Rate_parity_Set_M,snr_M]=meshgrid(Rate_parity_Set,SNRLin);
    FInterpolation2 = scatteredInterpolant(10.*log10(snr_M(:)),Rate_parity_Set_M(:),10.*log10(WEP_set(:)));
    WEP= 10.^(FInterpolation2(10.*log10(gamma),parity_rate)./10);
    return
    
elseif transmission_rate==3
    
    load WEP_Turbo_16QAM_rate3poincenage6_3gpp.mat
    SNRLin(SNRLin==0)=1e-50;
    WEP_set(isnan(WEP_set))=1e-50;
    WEP_set(WEP_set==0)=1e-50;
    [Rate_parity_Set_M,snr_M]=meshgrid(Rate_parity_Set,SNRLin);
    FInterpolation3 = scatteredInterpolant(10.*log10(snr_M(:)),Rate_parity_Set_M(:),10.*log10(WEP_set(:)));
    WEP= 10.^(FInterpolation3(10.*log10(gamma),parity_rate)./10);
    return

elseif transmission_rate==3.75
    
    load WEP_Turbo_16QAM_rate3.75poincenage6_3gpp.mat
    SNRLin(SNRLin==0)=1e-50;
    WEP_set(isnan(WEP_set))=1e-50;
    WEP_set(WEP_set==0)=1e-50;
    [Rate_parity_Set_M,snr_M]=meshgrid(Rate_parity_Set,SNRLin);
    FInterpolation4 = scatteredInterpolant(10.*log10(snr_M(:)),Rate_parity_Set_M(:),10.*log10(WEP_set(:)));
    WEP= 10.^(FInterpolation4(10.*log10(gamma),parity_rate)./10);
    return
 
else
    
    error( 'the .mat is not available for this transmission rate to estimate the WEP' );
    
end
    
    
    
    
    


