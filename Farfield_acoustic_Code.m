%% Title : Code for pressure data
%% Coded by: Kaurab Gautam
%% Input required
clc;
clear all;
close all; 
T_aq = 2 %sec
F_aq = 204800 % Acquizition frequency or sampling rate
N = T_aq*F_aq
n_mics = 16 % no of mics

fullname = fullfile('C:\Users\Kaurab Gautam\OneDrive - University of Cincinnati\Desktop\MyThesis\TestData\FarField\TR1p0\TwinMajor\Test10_19\WithShroud\NPR4p0', 'Pressure_FF_Pos001.dat');
data = load(fullname);
data(:,1) = [];
pressure = data;

%%
pre_fluc(n_mics,N) =0 ;
OASPL_rms(n_mics,1) = 0;
% ## Normalizing the value of P and assigning new vwctor Pp
sum = 0
for j = 1: length(n_mics)
    p_mean=mean(pressure(j,:)); % makes average of each rows ( which has 1024000 data vales)
end
pre_fluc =pressure -p_mean;
OASPL_rms = rms(pre_fluc,2);

% This gives OASPL in pa
%% Calculation OASPL in db
P_ref = 20e-6 % referecne pressure
OASPL_rms = 20*log10(OASPL_rms./P_ref)
angles = [152, 148, 144, 140, 136, 132, 126, 120, 116, 110, 105,100, 90, 70, 60, 45]

%% Plotting
figure;
plot(angles,OASPL_rms,"r", 'LineWidth',1.1)
xlim([45 155])
ylim([105 135])
grid on
xlim([45 155])
ylim([105 140])
grid on
txt = [' $ \psi $ degree']
title(txt, 'interpreter','latex','fontsize',11)
ylabel('OASPL, dB','interpreter','latex','fontsize',11)
xlabel(' Observation angles ( $ \psi $ ), Hz','interpreter','latex','fontsize',11)

%% Initial conditions for SPL data
blk_size = 4096  % this is block size which is defined as bs = N/nb where nb is number of block = 200
nb = N/blk_size;
F_res = F_aq/blk_size; % Frequency data resolution per second
ffi = F_res:F_res:F_aq; % integration domain
Fb(blk_size, (N*F_res/F_aq)) = 0 ;  %Fb = np.zeros(bs,nb) Coz we know nb = N*(df/fs))
SPL(blk_size/2,n_mics) = 0 ; %Calculating SPL for each frequency
OASPL_rms2(1,n_mics) = 0 ;

Fdom = ffi(1:1:(blk_size/2)) ; % Frequency domain
%% SPL calculation section
for u = 1:n_mics
    X = pre_fluc(u,:);
    b = 0
    for k = 1:(N*F_res)/(F_aq)
        a = b+1;
        b=k*blk_size;
        Fb(:,k)=abs(fft(X(a:b))./blk_size);
    end
    Fb = Fb.^2;
    F = mean(Fb,2);
    PSDB = (1/F_res)*F;
    q = trapz(ffi,PSDB);
    OASPL_2(1,u) = 10*log10(q/P_ref^2);
    SPL(:,u)=10*log10(((2*F_res)*PSDB(1:blk_size/2,1))/(P_ref^2));
end
%% Plotting SPLs
figure;
for i = 1:n_mics/2
    ax(i)= subplot(2,4,i)
    semilogx(Fdom,SPL(:,i))
    xlim([300 150000])
    ylim([60 125])
    grid on
    txt = [' $ \psi $ =', num2str(angles(i)), ' degree']
    title(txt, 'interpreter','latex','fontsize',11)
    ylabel('SPL, dB','interpreter','latex','fontsize',11)
    xlabel('Frequency, Hz','interpreter','latex','fontsize',11)
end
figure;
for i = 1:n_mics/2
    ax(i)= subplot(2,4,i)
    semilogx(Fdom,SPL(:,i+(n_mics/2)))
    xlim([300 150000])
    ylim([60 125])
    grid on
    txt = [' $ \psi $ =', num2str(angles(i+(n_mics/2))), ' degree']
    title(txt, 'interpreter','latex','fontsize',11)
    ylabel('SPL, dB','interpreter','latex','fontsize',11)
    xlabel('Frequency, Hz','interpreter','latex','fontsize',11)
end
%% Saving data
save('2SNPR4p0','SPL',"OASPL_rms")
save('Fdom','Fdom')

