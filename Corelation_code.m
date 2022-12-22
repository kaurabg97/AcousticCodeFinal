%% Title : Corelation code
%% Coded by: Kaurab Gautam
%% Input required
clc;
clear all;
close all;
T_aq = 10 % sec. This is the time for which sample is taken
F_aq = 204800 % Acquizition frequency of the mic
%N = T_aq*F_aq % Total data collected by each mics

fullname = fullfile('C:\Users\Kaurab Gautam\OneDrive - University of Cincinnati\Desktop\MyThesis\TestData\FarField\TR1p0\TwinMajor\Test10_27\NoShroud\NPR4p0', 'Pressure_FF_Pos001.dat');
%fullname = fullfile('C:\Users\Kaurab Gautam\OneDrive - University of Cincinnati\Desktop\MyThesis\TestData\FarField\TR1p0\TwinMajor\AatreshFF\Major\TR_1p0\NPR_4p0','Pressure_FF_Pos001.dat');


data = load(fullname);
data(:,1) = [];
pressure = data;
N_mics = size(pressure,1) % no of mics
N = size(pressure,2) % Total data collected by each mics
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
OASPL = 20*log10(Pressure_rms./P_ref);


%% Initial conditions for SPL data

blk_size = 4096 ;  % this is block size whicah is defined as bs = N/nb where nb is number of block = 200
F_res = F_aq/blk_size; % Frequency data resolution
ffi = F_res:F_res:F_aq; % integration domain % does it start with zero?
Fb(blk_size, (N*F_res/F_aq)) = 0; %Fb = (blk size , number of blocks)
y(blk_size, (N*F_res/F_aq)) = 0; % To store phase value %Calculating SPL for each frequency, blk size/2 coz it is symmetric
Fdom = ffi(1:1:(blk_size/2)); % Frequency domain
SPL(N_mics,blk_size) = 0;

%%
figure;
[acfNS,lagsNS]= xcorr(pressure(1,:), 'normalized');
plot(lagsNS/F_aq, acfNS)
xlim([-0.0005 0.0005])

%%
% [xc,lags] = autocorr(pre_fluc(1,:));
% figure;
% plot(lags/F_aq, xc)
% %%
% figure;
% [acf2,lags2]= crosscorr(pre_fluc(1,:), pre_fluc(2,:),NumLags=50,NumSTD=1.96);
% plot(lags2, acf2)


% % % % %%
% plot(lagsAA/F_aq, acfAA, 'r', LineWidth=1.5)
% hold on
% plot(lagsWS/F_aq, acfWS,'k',LineWidth=1.5)
% hold on
% plot(lagsNS/F_aq, acfNS, 'b',LineWidth=1.5)
% xlim([-0.0005 0.0005])
% xlabel('Time, [s]','interpreter','latex','fontsize',11)
% ylabel('Autocorelation','interpreter','latex','fontsize',11)
% title( 'Autocorelation result for NPR 4.0 (152 degree mic) ','interpreter','latex','fontsize',11)
% legend('Aatresh-4.0', 'Shroud-4.0', 'No shroud 4.0', 'interpreter','latex','fontsize',11)