%SPL Comparision
clearvars

%% OSU data
root = 'C:\Users\Kaurab Gautam\OneDrive - University of Cincinnati\Desktop\MyThesis\TestData\OSU\050272022_FF_AA90\';
OSUD = load([root,'OSU4p0.MAT']);
OSUDSPL = OSUD.SPL;
OSUDFdom = OSUD.Fdom;


%% My data

root = 'C:\Users\Kaurab Gautam\OneDrive - University of Cincinnati\Desktop\MyThesis\TestData\FarField\TR1p0\TwinMajor\Test10_27\NoShroud\';
MD = load([root,'BaselineNPR4p0.MAT']);
MDSPL = MD.SPL;
%MSFdom = MD.Fdom;

%% Aatresh data

root = 'C:\Users\Kaurab Gautam\OneDrive - University of Cincinnati\Desktop\MyThesis\TestData\ToCompare\';
AD = load([root,'ANPR4p0.MAT']);
ADFdom = load([root,'Fdom.MAT']);
MDFdom = ADFdom;
ADSPL = AD.SPL;
%ADFdom = AD.Fdom

%% plotitng
% for 150 degree
figure(1)
subplot(131)
semilogx(OSUDFdom(1,:),OSUDSPL(:,1),'red', LineWidth =1.3)
%ylim([60 125])
grid on
%title(txt, 'interpreter','latex','fontsize',9)
ylabel('SPL, [dB]','interpreter','latex','fontsize',9)
xlabel('Frequency, [HZ]','interpreter','latex','fontsize',9)
legend

hold on

semilogx(MDFdom.Fdom(1,:),MDSPL(:,1),'green', LineWidth=1.3)
%ylim([60 125])
grid on
%title(txt, 'interpreter','latex','fontsize',9)
ylabel('SPL, [dB]','interpreter','latex','fontsize',9)
xlabel('Frequency, [HZ]','interpreter','latex','fontsize',9)
legend
hold on

semilogx(ADFdom.Fdom(1,:),ADSPL(:,1),'black',LineWidth=1.3)
%ylim([60 125])
grid on
%title(txt, 'interpreter','latex','fontsize',9)
ylabel('SPL, [dB]','interpreter','latex','fontsize',9)
xlabel('Frequency, [HZ]','interpreter','latex','fontsize',9)
legend
title("SPL comparisions at 150 degree mic ",'interpreter','latex','fontsize',11)
xlim([300 104800])
ylim([60 130])
legend('OSU result','My baseline','Aatresh result','interpreter','latex','fontsize',11)
pbaspect([1 1 1])

%% for 120 degree
subplot(132)
semilogx(OSUDFdom(1,:),OSUDSPL(:,6),'red', LineWidth =1.3)
%ylim([60 125])
grid on
%title(txt, 'interpreter','latex','fontsize',9)
ylabel('SPL, [dB]','interpreter','latex','fontsize',9)
xlabel('Frequency, [HZ]','interpreter','latex','fontsize',9)
legend

hold on

semilogx(MDFdom.Fdom(1,:),MDSPL(:,8),'green', LineWidth=1.3)
%ylim([60 125])
grid on
%title(txt, 'interpreter','latex','fontsize',9)
ylabel('SPL, [dB]','interpreter','latex','fontsize',9)
xlabel('Frequency, [HZ]','interpreter','latex','fontsize',9)
legend
hold on

semilogx(ADFdom.Fdom(1,:),ADSPL(:,8),'black',LineWidth=1.3)
%ylim([60 125])
grid on
%title(txt, 'interpreter','latex','fontsize',9)
ylabel('SPL, [dB]','interpreter','latex','fontsize',9)
xlabel('Frequency, [HZ]','interpreter','latex','fontsize',9)
legend
title("SPL comparisions at 120 degree mic ",'interpreter','latex','fontsize',11)
xlim([300 104800])
ylim([60 130])
legend('OSU result','My baseline','Aatresh result','interpreter','latex','fontsize',11)
pbaspect([1 1 1])

%%
subplot(133)
semilogx(OSUDFdom(1,:),OSUDSPL(:,8),'red', LineWidth =1.3)
%ylim([60 125])
grid on
%title(txt, 'interpreter','latex','fontsize',9)
ylabel('SPL, [dB]','interpreter','latex','fontsize',9)
xlabel('Frequency, [HZ]','interpreter','latex','fontsize',9)
legend

hold on

semilogx(MDFdom.Fdom(1,:),MDSPL(:,13),'green', LineWidth=1.3)
%ylim([60 125])
grid on
%title(txt, 'interpreter','latex','fontsize',9)
ylabel('SPL, [dB]','interpreter','latex','fontsize',9)
xlabel('Frequency, [HZ]','interpreter','latex','fontsize',9)
legend
hold on

semilogx(ADFdom.Fdom(1,:),ADSPL(:,13),'black',LineWidth=1.3)
%ylim([60 125])
grid on
%title(txt, 'interpreter','latex','fontsize',9)
ylabel('SPL, [dB]','interpreter','latex','fontsize',9)
xlabel('Frequency, [HZ]','interpreter','latex','fontsize',9)
legend
title("SPL comparisions for 90 degree mic ",'interpreter','latex','fontsize',11)
xlim([300 104800])
ylim([60 130])
legend('OSU result','My baseline','Aatresh result','interpreter', 'latex','fontsize',11)
pbaspect([1 1 1])
