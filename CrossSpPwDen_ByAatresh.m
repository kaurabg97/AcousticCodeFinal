%% Title : Near Field Code
%% Coded by: Kaurab Gautam
%% Input required
clc; clearvars;
T_aq = 10 % sec. This is the time for which sample is taken
F_aq = 204800 % Acquizition frequency of the mic
Te = 20  % Degree C
De= 0.0184658 % m
NPR= 4.0;
gamma = 1.4;
T = 293.65 % [K]
R = 287.052874 % [J/Kg.K]
Ve = sqrt(2*gamma*R*T*((NPR)^((gamma-1)/gamma)-1)/(gamma-1)) % [m/s]

fullname = fullfile('C:\Users\Kaurab Gautam\OneDrive - University of Cincinnati\Desktop\MyThesis\TestData\NearField\Test11-13\NoShroudV3\NPR2p5', 'Pressure_FF_Pos001.dat');
data = load(fullname);
data(:,1) = [];
pressure = data;
%data= data(:,1024001:1433600);
N_mics = size(pressure,1) % no of mics
N = size(pressure,2)
pre_fluc(N_mics,N) =0 ; % pre_fluc=zeros(size(pressure))
Pressure_rms(N_mics,1) = 0;
% ## Normalizing the value of P and assigning new vwctor Pp
% makes average of each columns i.e for each microphoen it takes mean of all 2048100 value ( which has 1024000 data vales)
Pm=mean(pressure,2);
pre_fluc =pressure -Pm;
Pressure_rms = rms(pre_fluc,2)
%% Calculation OASPL in db
P_ref = 20e-6 % referecne pressure
OASPL = 20*log10(Pressure_rms./P_ref)

%% Initial conditions for SPL data

blk_size = 4096 ;  % this is block size whicah is defined as bs = N/nb where nb is number of block = 200
F_res = F_aq/blk_size; % Frequency data resolution
ffi = F_res:F_res:F_aq; % integration domain % does it start with zero?
Fb(blk_size, (N*F_res/F_aq)) = 0; %Fb = (blk size , number of blocks)
y(blk_size, (N*F_res/F_aq)) = 0; % To store phase value %Calculating SPL for each frequency, blk size/2 coz it is symmetric
Fdom = ffi(1:1:(blk_size/2)); % Frequency domain
SPL(N_mics,blk_size) = 0;

%% For shroud case
%% Input required
T_aq = 10 % sec. This is the time for which sample is taken
F_aq = 204800 % Acquizition frequency of the mic

fullname1 = fullfile('C:\Users\Kaurab Gautam\OneDrive - University of Cincinnati\Desktop\MyThesis\TestData\NearField\Test11-13\WithShroud\TR1p0\NPR2p5', 'Pressure_FF_Pos001.dat');
data1 = load(fullname1);
data1(:,1) = [];
pressure1 = data1;
%data1= data1(:,1024001:1433600);
N_mics = size(pressure1,1) % no of mics
N = size(pressure1,2)
pre_fluc1(N_mics,N) =0 ; % pre_fluc=zeros(size(pressure))
Pressure_rms1(N_mics,1) = 0;
% ## Normalizing the value of P and assigning new vwctor Pp
% makes average of each columns i.e for each microphoen it takes mean of all 2048100 value ( which has 1024000 data vales)
Pm1=mean(pressure1,2);
pre_fluc1 =pressure1 -Pm1;
Pressure_rms1 = rms(pre_fluc1,2)
%% Calculation OASPL in db
P_ref = 20e-6 % referecne pressure
OASPL1 = 20*log10(Pressure_rms1./P_ref)

% for 1 and 5
% namig is like , Fr of baseline between 1 and 5 = FrB1
[croSpecB1,FrB1] = cpsd(pre_fluc(3,:),pre_fluc(6,:),[],[],4096,204800);
[croSpecS1,FrS1] = cpsd(pre_fluc1(3,:),pre_fluc1(6,:),[],[],4096,204800);
phsDifB1 = angle(croSpecB1);
phsDifS1 = angle(croSpecS1);

% for mic 2 and 4
[croSpecB2,FrB2] = cpsd(pre_fluc(4,:),pre_fluc(5,:),[],[],4096,204800);
[croSpecS2,FrS2] = cpsd(pre_fluc1(4,:),pre_fluc1(5,:),[],[],4096,204800);
phsDifB2 = angle(croSpecB2);
phsDifS2 = angle(croSpecS2);
%% Code by Aatresh
C = input('Frequency or Strouhaoul number [1/2] ? ')
if C ==1
    figure(3);
    subplot(2,1,1)
    %power plot
    semilogx((FrB1+50),abs(croSpecB1),'LineWidth',1.2,'Color','r')
    hold on
    semilogx((FrS1+50),abs(croSpecS1),'LineWidth',1.2,'Color','k')
    title("Magnitude of CPSD between mics 3 and 6",'interpreter','latex','fontsize',11)
    xlabel('Frequency, [Hz]','interpreter','latex','fontsize',11)
    ylabel('Power, [W/Hz]','interpreter','latex','fontsize',11)
    legend('Baseline case','Shroud case','interpreter','latex','fontsize',11)
    xlim([4500, 25000])
    %ylim([0, 1350])
    %phase plot
    subplot(2,1,2)
    semilogx((FrB1+50),phsDifB1,'o','Color','r')
    hold on
    semilogx((FrS1+50),phsDifS1,'o','Color','k')
    title("Phase difference, [radians] ",'interpreter','latex','fontsize',11)
    xlabel('Frequency, [Hz]','interpreter','latex','fontsize',11)
    ylabel('Phase difefrence, Radian','interpreter','latex','fontsize',11)
    ylim([-pi pi])
    legend('Baseline case','Shroud case','interpreter','latex','fontsize',11)
    xlim([4500, 22000])
    %% for 2 and 4
    figure(4);
    subplot(2,1,1)
    %power plot
    semilogx((FrB2+50),abs(croSpecB2),'LineWidth',1.2,'Color','r')
    hold on
    semilogx((FrB2+50),abs(croSpecS2),'LineWidth',1.2,'Color','k')
    title("Magnitude of CPSD between mics 4 and 5",'interpreter','latex','fontsize',11)
    xlabel('Frequency, [Hz] ','interpreter','latex','fontsize',11)
    ylabel('Power, [W/Hz] ','interpreter','latex','fontsize',11)
    xlim([4500, 22000])
    legend('Baseline case','Shroud case','interpreter','latex','fontsize',11)
    %phase plot
    subplot(2,1,2)
    semilogx((FrB2+50),phsDifB2,'o','Color','r')
    hold on
    semilogx((FrB2+50),phsDifS2,'o','Color','k')
    title("Phase difference",'interpreter','latex','fontsize',11)
    xlabel('Frequency, [Hz]','interpreter','latex','fontsize',11)
    ylabel('Phase difefrence, [radians]','interpreter','latex','fontsize',11)
    ylim([-pi pi])
    xlim([4500, 22000])
    legend('Baseline case','Shroud case','interpreter','latex','fontsize',11)
else
    figure(3);
    %power plot
    St1 = (((FrB1+50)*De)/Ve);
    St2 = (((FrB2+50)*De)/Ve);

    subplot(2,1,1)
    semilogx(St1,abs(croSpecB1),'LineWidth',1.2,'Color','r')
    hold on
    logx(St2,abs(croSpecS1),'LineWidth',1.2,'Color','k')
    title("Magnitude of CPSD between mics 1 and 5",'interpreter','latex','fontsize',11)
    xlabel('St','interpreter','latex','fontsize',11)
    ylabel('Power, [W/Hz]','interpreter','latex','fontsize',11)
    legend('Baseline case','Shroud case','interpreter','latex','fontsize',11)
    xlim([0, 4.6])
    ylim([0, 360])
    %phase plot
    subplot(2,1,2)
    semilogx(St1,phsDifB1,'LineWidth',1.2,'Color','r')
    hold on
    semilogx(St2,phsDifS1,'LineWidth',1.2,'Color','k')
    title("Phase difference",'interpreter','latex','fontsize',11)
    xlabel('St','interpreter','latex','fontsize',11)
    ylabel('Phase difefrence, [radians]','interpreter','latex','fontsize',11)
    ylim([-pi pi])
    legend('Baseline case','Shroud case','interpreter','latex','fontsize',11)
    xlim([0, 4.6])
    %% for 2 and 4
    figure(4);
    %power plot
    subplot(2,1,1)
    semilogx(St1,abs(croSpecB2),'LineWidth',1.2,'Color','r')
    hold on
    semilogx(St2,abs(croSpecS2),'LineWidth',1.2,'Color','k')
    title("Magnitude of CPSD between mics 2 and 4",'interpreter','latex','fontsize',11)
    xlabel('St','interpreter','latex','fontsize',11)
    ylabel('Power, W/Hz','interpreter','latex','fontsize',11)
    xlim([0, 4.6])
    ylim([0, 360])
    legend('Baseline case','Shroud case','interpreter','latex','fontsize',11)
    %phase plot
    subplot(2,1,2)
    semilogx(St1,phsDifB2,'LineWidth',1.2, 'Color','r')
    hold on
    semilogx(St2,phsDifS2,'LineWidth',1.2,'Color','k')
    title("Phase difference",'interpreter','latex','fontsize',11)
    xlabel('St','interpreter','latex','fontsize',11)
    ylabel('Phase difefrence, Radian','interpreter','latex','fontsize',11)
    ylim([-pi pi])
    xlim([0, 4.6])
    legend('Baseline case','Shroud case','interpreter','latex','fontsize',11)

end

