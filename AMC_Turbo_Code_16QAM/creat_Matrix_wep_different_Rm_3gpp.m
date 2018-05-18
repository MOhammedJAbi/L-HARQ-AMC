clc;clear;

load perSbrqTurboLte.mat


rate=3.75
poincenage=3
ind_rate=find(Rset==rate)
max_ind_rate_parity=length(Rpset)-1
ind_poin=find(Rmset==poincenage)

% i=30;
% 
% length(snr{16,i}(1,:))

SNRLin(1)=0;
WEP_set(1,(1:length(Rpset)))=1;

for i=1:length(snr{ind_rate,1}(1,:))
    
    SNRLin(i+1)=snr{ind_rate,1}(1,i);
    
    for j=1:max_ind_rate_parity
        Rate_parity_Set(j)=Rpset(j)*rate/32;
        WEP_set(i+1,j)=result{ind_rate,j,ind_poin}.perCum(2,i)./result{ind_rate,j,ind_poin}.experiments(i);
    end
end

SNRLin(i+2)=1e100;
WEP_set(i+2,(1:ind_rate))=0;
Rate_parity_Set(j+1)=Rpset(j+1)*rate/32;
WEP_set(:,j+1)=0;

save(['WEP_Turbo_16QAM_rate',num2str(rate),'poincenage6_3gpp.mat'],'SNRLin','Rate_parity_Set','WEP_set')