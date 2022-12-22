%% Title : Near Field Code
%% Coded by: Kaurab Gautam
%% Input required
clc;
clear all;
close all;
T_aq = 10 % sec. This is the time for which sample is taken
F_aq = 204800 % Acquizition frequency of the mic
N = T_aq*F_aq % Total data collected by each mics

fullname = fullfile('C:\Users\Kaurab Gautam\OneDrive - University of Cincinnati\Desktop\MyThesis\TestData\NearField\MajorAxis\NPR2p5', 'Pressure_FF_Pos001.dat');
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
%angles = [152, 148, 144, 140, 136, 132, 126, 120, 116, 110, 105,100, 90, 70, 60, 45]


%% Initial conditions for SPL data

blk_size = 4096 ;  % this is block size whicah is defined as bs = N/nb where nb is number of block = 200
F_res = F_aq/blk_size; % Frequency data resolution
ffi = F_res:F_res:F_aq; % integration domain % does it start with zero?
Fb(blk_size, (N*F_res/F_aq)) = 0; %Fb = (blk size , number of blocks)
y(blk_size, (N*F_res/F_aq)) = 0; % To store phase value %Calculating SPL for each frequency, blk size/2 coz it is symmetric
Fdom = ffi(1:1:(blk_size/2)); % Frequency domain
SPL(N_mics,blk_size) = 0;
%% SPL calculation section
for u = 1:n_mics
    X = pre_fluc(u,:)
    b = 0
    for k = 1:(N*F_res)/(F_aq)
        a = b+1;
        b=k*blk_size;
        Fb(:,k)=fft(X(a:b)); % calculates FFT for each block
        Fba(:,k)=abs(Fb./blk_size); %v calculates absolute value
        y(:,k) = angle(Fb(:,k));
    end
    Fba = Fba.^2; % calculates square, dim = [ blk_size , nb)
    F = mean(Fba,2); % calculates row wise mean of the absoulte value of fourier transform dim = [ blk_size , 1]
    PSDB = (1/F_res)*F; % Normalizes F
    q = trapz(ffi,PSDB); % integration for PSDB ( i.e fourier amplitude) to the interval of ffi
    SPL(:,u)=10*log10(((2*F_res)*PSDB(1:blk_size/2,1))/(P_ref^2)); % calculates ffi for half of the block size
end
% %% \Plotting Phases
%% Phase shift through averaging
for mic = 1:N_mics
    X = pre_fluc(mic,:);
    b = 0
    for k = 1:(N*F_res)/(F_aq) % Loop goes from 1 to total number of block
        a = b+1;
        b=k*blk_size;
        Fbr(:,k)=fft(X(a:b)); % real fft 
        Fb(:,k)=abs(Fbr(:,k))/blk_size;
    end
    Fb = Fb.^2;
    F = mean(Fb,2);
    MeanFbr = mean(Fbr,2);
    PSDB = (1/F_res)*F;
    q = trapz(ffi,PSDB);
    OASPL_2(1,mic) = 10*log10(q/P_ref^2);
    SPL(:,mic)=10*log10(((2*F_res)*PSDB(1:blk_size/2,1))/(P_ref^2));
    FFTR(mic,:)=MeanFbr;
end

%% Phase shift without averaging
Wave1 = fft(pre_fluc(1,:))
Wave2 = fft(pre_fluc(2,:))
[value1,index1]=max(2*abs(Wave1(1:F_aq/2+1)));
pV1 = angle(Wave1(index1));
[value2,index2]=max(2*abs(Wave2(1:F_aq/2+1)));
pV2 = angle(Wave2(index2))
pd = pV1-pV2 
%% Plotting SPLs
figure;
angles = [1 2]
for i = 1:N_mics
    ax(i)= subplot(1,2,i)
    semilogx(Fdom,SPL(:,i))
    xlim([300 150000])
    ylim([60 125])
    grid on
    txt = [' Mic number =', num2str(angles(i))']
    title(txt, 'interpreter','latex','fontsize',11)
    ylabel('SPL, dB','interpreter','latex','fontsize',11)
    xlabel('Frequency, Hz','interpreter','latex','fontsize',11)
end

%%
mic1 = FFTR(1,:)
mic2 = FFTR(2,:)

[value,index]=max(2*abs(SPL((1:Fdom/2+1))));
pV1 = angle(mic1(index));
[value,index]=max(2*abs(SPL((2:Fdom/2+1))));
pV2 = angle(mic2(index));
% [value,index]=max(2*abs(FY2(1:NFFT/2+1)));
% pV2 = angle(FY2(index))
pd = pV1-pV2 
Ang =angle(FFTR)

%signal windowing



