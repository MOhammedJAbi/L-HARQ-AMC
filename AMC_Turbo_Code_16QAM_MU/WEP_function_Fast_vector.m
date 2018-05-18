function [ WEP] = WEP_function_Fast_vector(transmission_rate,parity_rate,gamma,FInterpolation1,FInterpolation2,FInterpolation3,FInterpolation4)

% This function estimates the WEP
if transmission_rate==1.5
    WEP= 10.^(FInterpolation1(10.*log10(gamma),parity_rate.*ones(size(gamma)))./10);
    WEP=min(max(WEP,0),1);
    return

elseif transmission_rate==2.25
    WEP= 10.^(FInterpolation2(10.*log10(gamma),parity_rate.*ones(size(gamma)))./10);
    WEP=min(max(WEP,0),1);
    return
    
elseif transmission_rate==3
    WEP= 10.^(FInterpolation3(10.*log10(gamma),parity_rate.*ones(size(gamma)))./10);
    WEP=min(max(WEP,0),1);
    return

elseif transmission_rate==3.75
    WEP= 10.^(FInterpolation4(10.*log10(gamma),parity_rate.*ones(size(gamma)))./10);
    WEP=min(max(WEP,0),1);
    return
 
else
    
    error( 'the .mat is not available for this transmission rate to estimate the WEP' );
    
end
    
    
    
    
    


