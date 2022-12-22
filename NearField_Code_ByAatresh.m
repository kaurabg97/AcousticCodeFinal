%% Title : Near Field Code
%% Coded by: Kaurab Gautam
%% Input required
clc;
clear all;
close all;
T_aq = 10 % sec. This is the time for which sample is taken
F_aq = 204800 % Acquizition frequency of the mic
N = T_aq*F_aq % Total data collected by each mics

fullname = fullfile('C:\Users\Kaurab Gautam\OneDrive - University of Cincinnati\Desktop\MyThesis\TestData\NearField\NearField10_28\TR1p0\V1\NPR3p0', 'Pressure_FF_Pos001.dat');
data = load(fullname);
data(:,1) = [];
pressure = data;
N_mics = size(pressure,1) % no of mics

%%
pre_fluc(N_mics,N) =0 ; % pre_fluc=zeros(size(pressure))
Pressure_rms(N_mics,1) = 0;
% ## Normalizing the value of P and assigning new vwctor Pp
% makes average of each columns i.e for each microphoen it takes mean of all 2048100 value ( which has 1024000 data vales)
Pm=mean(pressure,2);
pre_fluc =pressure -Pm;
Pressure_rms = rms(pre_fluc,2)

% This gives OASPL in pa
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

%% Code by Aatresh
[croSpec,Fr] = cpsd(pre_fluc(1,:),pre_fluc(2,:),[],[],4096,204800);
phsDif = angle(croSpec);
phsDif = angle(croSpec);
figure(1);
subplot(211)
semilogx(Fr,abs(croSpec),'LineWidth',1.2)
hold on
title("Magnitude of CS analysis at NPR 2.5",'interpreter','latex','fontsize',11)
xlabel('Frequency, Hz','interpreter','latex','fontsize',11)
ylabel('Power, W/Hz','interpreter','latex','fontsize',11)
subplot(212)
semilogx(Fr,phsDif,'LineWidth',1.2)
hold on
title("Phase difference at NPR 2.5",'interpreter','latex','fontsize',11)
xlabel('Frequency, Hz','interpreter','latex','fontsize',11)
ylabel('Phase difefrence, Radian','interpreter','latex','fontsize',11)
ylim([-pi pi])
%xlim([10000 11000])
%%
[s,f,t] = spectrogram(pre_fluc(1,:),4096,[],[],204800);
pcolor(t, f, abs(s).^2)
