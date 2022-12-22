%% Title : Near Field Code
%% Coded by: Kaurab Gautam
%% Input required
clc;
clear all;
close all;
T_aq = 2 % sec. This is the time for which sample is taken
F_aq = 204800 % Acquizition frequency of the mic
N = T_aq*F_aq % Total data collected by each mics

root = 'C:\Users\Kaurab Gautam\OneDrive - University of Cincinnati\Desktop\MyThesis\TestData\OSU\050272022_FF_AA90';

datafile = load([root,'\M160_pressure_data.mat'])
%data1 = datafile.M130_pressure.mat');
data = datafile.M160_pressure;
%% 
pressure = transpose(data);
N_mics = size(pressure,1) % no of mics

%%
pre_fluc=zeros(size(pressure));
Pressure_rms=zeros(N_mics,1);
% ## Normalizing the value of P and assigning new vwctor Pp
% makes average of each columns i.e for each microphoen it takes mean of all 2048100 value ( which has 1024000 data vales)
Pm=mean(pressure,2);
for i = 1:N_mics
 pre_fluc(i,:) =pressure(i,:) -Pm(i,1);
end
Pressure_rms = rms(pre_fluc,2);
% This gives OASPL in pa
%% Calculation OASPL in db
P_ref = 20e-6 % referecne pressure
OASPL = 20*log10(Pressure_rms./P_ref)
angles1 = [30, 35, 40, 45, 50, 60, 75,90];
angles = [150, 145, 140, 135, 130, 120, 105, 90]
%% Plotting
figure;
plot(angles,OASPL,"r", 'LineWidth',1.1)
xlim([45 155])
ylim([105 135])
grid on
xlim([90 150])
ylim([105 140])
grid on
txt = [' $ \psi $ degree']
title('OASPL-OSU at M = 1.60', 'interpreter','latex','fontsize',9)
ylabel('OASPL, dB','interpreter','latex','fontsize',7)
xlabel(' Observation angles ( $ \psi $ ), Hz','interpreter','latex','fontsize',7)
filename = 'OASPL_at_M1p60'; 
saveas(gca, fullfile(root, filename), 'png');
%% Initial conditions for SPL data

blk_size = 4096 ;  % this is block size whicah is defined as bs = N/nb where nb is number of block = 200
F_res = F_aq/blk_size; % Frequency data resolution
ffi = F_res:F_res:F_aq; % integration domain % does it start with zero?
Fb(blk_size, (N*F_res/F_aq)) = 0; %Fb = (blk size , number of blocks)
y(blk_size, (N*F_res/F_aq)) = 0; % To store phase value %Calculating SPL for each frequency, blk size/2 coz it is symmetric
Fdom = ffi(1:1:(blk_size/2)); % Frequency domain
SPL(blk_size/2,N_mics) = 0 ;
%% SPL calculation section
for u = 1:N_mics
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
for i = 1:N_mics
    ax(i)= subplot(2,4,i)
    semilogx(Fdom,SPL(:,i))
    xlim([300 150000])
    ylim([60 125])
    grid on
    txt = [' $ \psi $ =', num2str(angles(i)), ' degree']
    title(txt, 'interpreter','latex','fontsize',9)
    ylabel('SPL, dB','interpreter','latex','fontsize',7)
    xlabel('Frequency, Hz','interpreter','latex','fontsize',7)
end
filename = 'SPL_at_M1p60'; 
saveas(gca, fullfile(root, filename), 'png');
%% Saving data
matfile = fullfile(root, 'OSUAA90A1p60');
save(matfile,'SPL',"OASPL");

% save('OSUAA90A1p30','SPL',"OASPL_rms")
% % save('Fdom','Fdom')
## To change the modified status



