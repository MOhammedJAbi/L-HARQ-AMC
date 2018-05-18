clc;clear;clf;

figure(1)
semilogy(10,10,'-or','markerfacecolor','r','MarkerSize',7)
hold on
semilogy(10,10,'-sb','markerfacecolor','b','MarkerSize',7)
semilogy(10,10,'-^m','markerfacecolor','m','MarkerSize',7)
semilogy(10,10,'-vk','markerfacecolor','k','MarkerSize',7)

dop=0.05;
rho=(besselj(0,2*pi*dop)).^2;
lint = 4e2;
[x0, w] = GaussLegendre(lint);
pdf_joint=@(x,y,av_snr) exp(-(rho.*x+y)./((1-rho).*av_snr)).*besseli(0,2.*sqrt(rho.*x.*y)./((1-rho).*av_snr))./((1-rho).*av_snr);
av_snr=10.^(15./10);
load MI_16QAM.mat;
MI_16QAM=[0 MI_16QAM(1:end-1)' 4];
th=[0 10.^((SNR(1:end-1))'./10) 1e10];
Qam_th=@(R) interp1(MI_16QAM,th,R);

%%
load WEP_Turbo_16QAM_rate1.5poincenage6_3gpp.mat
WEP_set(isnan(WEP_set))=0;
WEP_set(WEP_set==0)=1e-100;
SNRLin(SNRLin==0)=1e-10;
SNRLindb=linspace(-10,25,200);
SNRLin_set=10.^(SNRLindb./10);

rate=1.5;
snr_th=Qam_th(rate);

for j=1:length(SNRLin_set)
    Dx=-log(10^-30).*(1-rho).*av_snr+rho.*SNRLin_set(j);
    xx = (x0+1)*Dx/2;
    yy=interp1(10.*log10(SNRLin),10.*log10(WEP_set(:,1)),10.*log10(xx));
    yy=10.^(yy./10);
    val=yy.*pdf_joint(SNRLin_set(j),xx,av_snr);
    WEP_set_tc(j)=sum(val.*w)*Dx/2;
    yy=xx<=Qam_th(rate);
    val=yy.*pdf_joint(SNRLin_set(j),xx,av_snr);
    WEP_set_p(j)=sum(val.*w)*Dx/2;
end

hold on
semilogy(10.*log10(SNRLin_set),WEP_set_tc,'-r','linewidth',1.2)

semilogy(10.*log10(SNRLin_set),WEP_set_p,'--r','linewidth',1.2)

p0=76;
p1=10;
y1=WEP_set_tc((p0:p1:end));
y11=WEP_set_p((p0:p1:end));
x1=SNRLin_set(p0:p1:end);
semilogy(10.*log10(x1),y1,'or','markerfacecolor','r','MarkerSize',7)
semilogy(10.*log10(x1),y11,'or','markerfacecolor','r','MarkerSize',7)

%%
load WEP_Turbo_16QAM_rate2.25poincenage6_3gpp.mat

WEP_set(isnan(WEP_set))=0;
WEP_set(WEP_set==0)=1e-100;
SNRLin(SNRLin==0)=1e-10;
SNRLindb=linspace(-10,25,200);
SNRLin_set=10.^(SNRLindb./10);

rate=2.25;
snr_th=Qam_th(rate);

for j=1:length(SNRLin_set)
    Dx=-log(10^-30).*(1-rho).*av_snr+rho.*SNRLin_set(j);
    xx = (x0+1)*Dx/2;
    yy=interp1(10.*log10(SNRLin),10.*log10(WEP_set(:,1)),10.*log10(xx));
    yy=10.^(yy./10);
    val=yy.*pdf_joint(SNRLin_set(j),xx,av_snr);
    WEP_set_tc(j)=sum(val.*w)*Dx/2;
    yy=xx<=Qam_th(rate);
    val=yy.*pdf_joint(SNRLin_set(j),xx,av_snr);
    WEP_set_p(j)=sum(val.*w)*Dx/2;
end

semilogy(10.*log10(SNRLin_set),WEP_set_tc,'-b','linewidth',1.2)

semilogy(10.*log10(SNRLin_set),WEP_set_p,'--b','linewidth',1.2)

p0=81;
p1=20;
y2=WEP_set_tc((p0:p1:end));
y22=WEP_set_p((p0:p1:end));
x2=SNRLin_set(p0:p1:end);
semilogy(10.*log10(x2),y2,'sb','markerfacecolor','b','MarkerSize',7)
semilogy(10.*log10(x2),y22,'sb','markerfacecolor','b','MarkerSize',7)
%%

load WEP_Turbo_16QAM_rate3poincenage6_3gpp.mat

WEP_set(isnan(WEP_set))=0;
WEP_set(WEP_set==0)=1e-100;
SNRLin(SNRLin==0)=1e-10;
SNRLindb=linspace(-10,25,200);
SNRLin_set=10.^(SNRLindb./10);

rate=3;
snr_th=Qam_th(rate);

for j=1:length(SNRLin_set)
    Dx=-log(10^-30).*(1-rho).*av_snr+rho.*SNRLin_set(j);
    xx = (x0+1)*Dx/2;
    yy=interp1(10.*log10(SNRLin),10.*log10(WEP_set(:,1)),10.*log10(xx));
    yy=10.^(yy./10);
    val=yy.*pdf_joint(SNRLin_set(j),xx,av_snr);
    WEP_set_tc(j)=sum(val.*w)*Dx/2;
    yy=xx<=Qam_th(rate);
    val=yy.*pdf_joint(SNRLin_set(j),xx,av_snr);
    WEP_set_p(j)=sum(val.*w)*Dx/2;
end

semilogy(10.*log10(SNRLin_set),WEP_set_tc,'-m','linewidth',1.2)

semilogy(10.*log10(SNRLin_set),WEP_set_p,'--m','linewidth',1.2)

p0=105;
p1=20;
y2=WEP_set_tc((p0:p1:end));
y22=WEP_set_p((p0:p1:end));
x2=SNRLin_set(p0:p1:end);
semilogy(10.*log10(x2),y2,'^m','markerfacecolor','m','MarkerSize',7)
semilogy(10.*log10(x2),y22,'^m','markerfacecolor','m','MarkerSize',7)
%%
load WEP_Turbo_16QAM_rate3.75poincenage6_3gpp.mat

WEP_set(isnan(WEP_set))=0;
WEP_set(WEP_set==0)=1e-100;
SNRLin(SNRLin==0)=1e-10;
SNRLindb=linspace(-10,25,200);
SNRLin_set=10.^(SNRLindb./10);

rate=3.75;
snr_th=Qam_th(rate);

for j=1:length(SNRLin_set)
    Dx=-log(10^-30).*(1-rho).*av_snr+rho.*SNRLin_set(j);
    xx = (x0+1)*Dx/2;
    yy=interp1(10.*log10(SNRLin),10.*log10(WEP_set(:,1)),10.*log10(xx));
    yy=10.^(yy./10);
    val=yy.*pdf_joint(SNRLin_set(j),xx,av_snr);
    WEP_set_tc(j)=sum(val.*w)*Dx/2;
    yy=xx<=Qam_th(rate);
    val=yy.*pdf_joint(SNRLin_set(j),xx,av_snr);
    WEP_set_p(j)=sum(val.*w)*Dx/2;
end

semilogy(10.*log10(SNRLin_set),WEP_set_tc,'-k','linewidth',1.2)

semilogy(10.*log10(SNRLin_set),WEP_set_p,'--k','linewidth',1.2)

p0=127;
p1=20;
y2=WEP_set_tc((p0:p1:end));
y22=WEP_set_p((p0:p1:end));
x2=SNRLin_set(p0:p1:end);
semilogy(10.*log10(x2),y2,'vk','markerfacecolor','k','MarkerSize',7)
semilogy(10.*log10(x2),y22,'vk','markerfacecolor','k','MarkerSize',7)
%%
xlabel('SNR')
ylabel('WEP')

grid

legend('R=1.5XX','R=2.25','R=3','R=3.75','Location','SouthWest')
axis([02 20 5e-4 1])
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




name_figure = 'Wep_TR';
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

