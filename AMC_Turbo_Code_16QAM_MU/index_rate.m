function [index] = index_rate(transmission_rate)

% This function returns the index of the transmission rate

if transmission_rate==1.5
    
    index=1;
    return

elseif transmission_rate==2.25
    
    index=2;
    return
    
elseif transmission_rate==3
    
    index=3;
    return

elseif transmission_rate==3.75
    
    index=4;
    return
 
else
    
    error( 'this rate is not available' );
    
end