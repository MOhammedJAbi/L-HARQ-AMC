
function out=WEP(snr,snr_th,g)
% if snr_th==0
%     out=0;
%     return
% end
out=exp(-g.*(snr./snr_th-1));
out(find(snr-snr_th<0))=1;   
out(find(snr_th==0))=0; 
end