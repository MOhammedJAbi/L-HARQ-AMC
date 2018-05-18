clc;clear;clf;

semilogy(10,10,'-^m','markerfacecolor','r','MarkerSize',7)
hold on
semilogy(10,10,'-sb','markerfacecolor','b','MarkerSize',7)
semilogy(10,10,'-or','markerfacecolor','m','MarkerSize',7)
semilogy(10,10,'-vk','markerfacecolor','k','MarkerSize',7)
semilogy(10,10,'-pg','markerfacecolor','g','MarkerSize',7)

%%
load WEP_Turbo_16QAM_rate3poincenage6_3gpp.mat

i=1
Rate_parity_Set(i)

semilogy(10.*log10(SNRLin),WEP_set(:,i),'-m','linewidth',1.2)

hold on

p0=104;
p1=2;
y1=WEP_set((p0:p1:end),1);
x1=SNRLin(p0:p1:end);
semilogy(10.*log10(x1),y1,'^m','markerfacecolor','m','MarkerSize',7)

%%
i=4
Rate_parity_Set(i)

semilogy(10.*log10(SNRLin),WEP_set(:,i),'-b','linewidth',1.2)
p0=96;
p1=2;
y2=WEP_set((p0:p1:end),i);
x2=SNRLin(p0:p1:end);
semilogy(10.*log10(x2),y2,'sb','markerfacecolor','b','MarkerSize',7)

%%

i=8
Rate_parity_Set(i)

semilogy(10.*log10(SNRLin),WEP_set(:,i),'-r','linewidth',1.2)
p0=86;
p1=3;
y2=WEP_set((p0:p1:end),i);
x2=SNRLin(p0:p1:end);
semilogy(10.*log10(x2),y2,'or','markerfacecolor','r','MarkerSize',7)

%%

i=11
Rate_parity_Set(i)

semilogy(10.*log10(SNRLin),WEP_set(:,i),'-k','linewidth',1.2)
p0=76;
p1=4;
y2=WEP_set((p0:p1:end),i);
x2=SNRLin(p0:p1:end);
semilogy(10.*log10(x2),y2,'vk','markerfacecolor','k','MarkerSize',7)

%%
i=15
Rate_parity_Set(i)
semilogy(10.*log10(SNRLin),WEP_set(:,i),'-g','linewidth',1.2)
p0=46;
p1=4;
y2=WEP_set((p0:p1:end),i);
x2=SNRLin(p0:p1:end);
semilogy(10.*log10(x2),y2,'pg','markerfacecolor','g','MarkerSize',7)

%%
xlabel('SNR')
ylabel('WEP')

grid

legend('R=0.0XX','R=0.5','R=1.3','R=1.8','R=2.6','Location','SouthWest')
axis([0 11.5 5e-4 1])
set(gca,'fontsize',11)
ll=legend;
set(ll,'FontSize',15);
%set(gcf,'PaperPositionMode','manual');

set(gcf,'PaperUnits','centimeters');
lar=512;
lon=2*lar/(1.1+sqrt(5));
set(gcf,'Position',[384 874 lar lon]);
set(gcf,'PaperPosition',[1.19 16 1.19 16]);
set(gca,'fontName','Times')
set(gcf,'PaperPositionMode','auto');




name_figure = 'Wep_BR';
print('-depsc','-r300',name_figure)



%load Turbo2_BRQ6.mat

% figure(1)
% 
% load WEP_Turbo_16QAM_rate1.5poincenage6_3gpp.mat
% 
% i=1;
% 
% Rate_parity_Set(i)
% semilogy(10.*log10(SNRLin),WEP_set(:,i),'-r')
% 
% hold on
% 
% i=5
% Rate_parity_Set(i)
% loglog(SNRLin,WEP_set(:,i),'-.m')
% 
% 
% i=10
% Rate_parity_Set(i)
% loglog(SNRLin,WEP_set(:,i),'-.g')
% 
% i=15
% Rate_parity_Set(i)
% loglog(SNRLin,WEP_set(:,i),'-.b')
% 
% axis([-5 6 1e-4 1])
% 
% legend('Rp=0','Rp=0.37','Rp=0.84','Rp=1.3')
% 
% xlabel('SNR')
% ylabel('WEP')
% 
% grid

% 
% figure(2)
% 
% load WEP_Turbo_16QAM_rate1.5poincenage6_3gpp.mat
% 
% for i=1:length(Rate_parity_Set)-1
%     semilogy(10.*log10(SNRLin),WEP_set(:,i))
%     hold on
% end
% axis([-5 6 1e-4 1])
% title('R=1.5')
% 
% 
% figure(3)
% 
% load WEP_Turbo_16QAM_rate2.25poincenage6_3gpp.mat
% 
% for i=1:length(Rate_parity_Set)-1
%     semilogy(10.*log10(SNRLin),WEP_set(:,i))
%     hold on
% end
% axis([-5 10 1e-4 1])
% title('R=2.25')
% 
% figure(4)
% 
% load WEP_Turbo_16QAM_rate3poincenage6_3gpp.mat
% 
% for i=1:length(Rate_parity_Set)-1
%     semilogy(10.*log10(SNRLin),WEP_set(:,i))
%     hold on
% end
% axis([-5 15 1e-4 1])
% title('R=3')
% 
% figure(5)
% 
% load WEP_Turbo_16QAM_rate3.75poincenage6_3gpp.mat
% 
% for i=1:length(Rate_parity_Set)-1
%     semilogy(10.*log10(SNRLin),WEP_set(:,i))
%     hold on
% end
% axis([-5 20 1e-4 1])
% title('R=4')
% 
% figure(6)
% 
% load MI_16QAM.mat
% MI_16QAM=[0 MI_16QAM(1:end-1)' 4];
% th=[0 10.^((SNR(1:end-1))'./10) 1e10];
% Qam_th=@(R) interp1(MI_16QAM,th,R);
% Qam_r=@(gamma) interp1(th,MI_16QAM,gamma);
% 
% load WEP_Turbo_16QAM_rate1.5poincenage6_3gpp.mat
% 
% load WEP_Turbo_16QAM_rate1.5poincenage6_3gpp.mat
% %[val,nan_v]=max(isnan(WEP_set(:,1)));
% %[val,nan_v]=min((WEP_set(:,1))>=1e-4);
% %j=nan_v;
% %j=(nan_v-1)*(nan_v>1)+length(SNRLin)*(nan_v==1);
% % SNRLin_m{i}=[0 SNRLin(1:j) 1e100];
% % WEP_m{i}=[1 (WEP_set((1:j),1))' 0];
% i=1;
% [C,IA,IC] = unique(WEP_set(:,1), 'last') ;
% output_wep = WEP_set((sort(IA)),1);
% output_SNRLin=SNRLin(sort(IA));
% [val,nan_v]=max(isnan(output_wep));
% j=(nan_v-1)*(nan_v>1)+length(IA)*(nan_v==1);
% SNRLin_m{i}=[1e-50 output_SNRLin(1:j) 1e100];
% WEP_m{i}=[1 output_wep(1:j)' 1e-50];
% 
% epsilon=1e-1;
% for i=1:length(Rate_parity_Set)-1
%     snr_th(i)=interp1(squeeze(WEP_m{1}(:,i)),10.*log10(squeeze(SNRLin_m{1})),epsilon);
% end
% plot(Rate_parity_Set(1:end-1),snr_th,'r',Rate_parity_Set(1:end-1),10.*log10(Qam_th(Rate_parity_Set(1:end-1))),'r')
% axis([-5 6 1e-4 1])
% title('R=1.5')
%     

