%% Title : Near Field Code
%% Coded by: Kaurab Gautam
%% Input required
clc; clearvars;
F_aq = 204800 % Acquizition frequency of the mic

%% Jet operating condition
Te = 20  % Degree C
De= 0.0184658 % m
NPR= 2.5;
gamma = 1.4;
T = 293.65 % [K]
R = 287.052874 % [J/Kg.K]
Ve = sqrt(2*gamma*R*T*((NPR)^((gamma-1)/gamma)-1)/(gamma-1)) % [m/s]
%% Data acquizition
%fullname = fullfile('C:\Users\Kaurab Gautam\OneDrive - University of Cincinnati\Desktop\MyThesis\TestData\FarField\TR1p0\TwinMajor\AatreshFF\Major\TR_1p0\NPR_2p5', 'Pressure_FF_Pos001.dat');
root = 'C:\Users\Kaurab Gautam\OneDrive - University of Cincinnati\Desktop\MyThesis\TestData\OSU\050272022_FF_AA90\';
data = load([root,'M155_pressure_data.MAT']);
data = data.M155_pressure;
%data(:,1) = [];
pressure = data;
pressure = transpose(pressure);
N_mics = size(pressure,1); % no of mics
N = size(pressure,2);
T_aq = N/F_aq;
%% Calculation OASPL in db
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
P_ref = 20e-6 % referecne pressure
OASPL = 20*log10(Pressure_rms./P_ref)
angles =[ 1, 2, 3, 4, 5, 6,7,8]
%angles = [152, 148, 144, 140, 136, 132, 126, 120, 116, 110, 105,100, 90, 70, 60, 45]

%% Plotting of OASPL
% figure;
% plot(angles,OASPL,"r", 'LineWidth',1.1)
% xlim([45 155])
% ylim([105 135])
% grid on
% xlim([45 155])
% ylim([105 140])
% grid on
% txt = [' $ \psi $ degree']
% title('OASPL at NPR 4.3 -shroud case', 'interpreter','latex','fontsize',11)
% ylabel('OASPL, dB','interpreter','latex','fontsize',11)
% xlabel(' Observation angles ( $ \psi $ ), Hz','interpreter','latex','fontsize',11)
% filename = 'OASPL_at_NPR4p3'; 
% saveas(gca, fullfile(root, filename), 'png');

%% Initial conditions for SPL data
blk_size = 4096 ;  % this is block size whicah is defined as bs = N/nb where nb is number of block = 200
F_res = F_aq/blk_size; % Frequency data resolution
ffi = F_res:F_res:F_aq; % integration domain % does it start with zero?
Fb(blk_size, (N*F_res/F_aq)) = 0; %Fb = (blk size , number of blocks)
y(blk_size, (N*F_res/F_aq)) = 0; % To store phase value %Calculating SPL for each frequency, blk size/2 coz it is symmetric
Fdom = ffi(1:1:(blk_size/2)); % Frequency domain
SPL(blk_size/2,N_mics) = 0 ;
% win = hann(N);
% win = win/sum(abs(win));
% 
% for i = N_mics
%     for j = N
%         pre_fluc(i,j) = pre_fluc(i,j).*win(j);
%     end
% end 
%% SPL calculation section
M = [2.72, 2.82, 3.01, 3.32, 3.58, 3.78, 3.35, 3.16]
for u = 1:N_mics
    X = pre_fluc(u,:);
    b = 0;
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
% correction for distance
    T = 1.135634 % 61.5 De or 44.7 in % target distance [m]
    M = [2.72, 2.82, 3.01, 3.32, 3.58, 3.78, 3.35, 3.16] % Measured distace [m]
    SPL(:,u) =  SPL(:,u)-20*log10(T/M(1,u));
end
% Plotting SPL
C = input('St or Hz? [1/2] : ')
if C == 1
    figure(2);
    for i = 1:N_mics
        %ax(i)= subplot(2,3,i)
        semilogx(((Fdom*De)/Ve),SPL(:,i),'linewidth',1.2,'DisplayName',num2str(i))
        xlim([0 4.6])
        %ylim([60 125])
        grid on
        txt = [' $ \psi $ =', num2str(angles(i)), ' degree']
        title(txt, 'interpreter','latex','fontsize',9)
        ylabel('SPL, [dB]','interpreter','latex','fontsize',9)
        xlabel('St','interpreter','latex','fontsize',9)
        hold on
        legend
    end
else
        figure(1);
    for i = 1:N_mics
        %ax(i)= subplot(2,3,i)
        semilogx(Fdom,SPL(:,i),'linewidth',1.2,'DisplayName',num2str(i))
        xlim([300 104800])
        %ylim([60 125])
        grid on
        txt = [' $ \psi $ =', num2str(angles(i)), ' degree']
        title(txt, 'interpreter','latex','fontsize',9)
        ylabel('SPL, [dB]','interpreter','latex','fontsize',9)
        xlabel('St','interpreter','latex','fontsize',9)
        hold on
        legend
    end
end

% filename = 'SPL_at_NPR4p3_part1'; 
% saveas(gca, fullfile(root, filename), 'png');
% figure;
% for i = 1:N_mics/2
%     ax(i)= subplot(2,4,i)
%     semilogx(Fdom,SPL(:,i+(N_mics/2)))
%     xlim([300 150000])
%     ylim([60 125])
%     grid on
%     txt = [' $ \psi $ =', num2str(angles(i+(N_mics/2))), ' degree']
%     title(txt, 'interpreter','latex','fontsize',9)
%     ylabel('SPL, dB','interpreter','latex','fontsize',9)
%     xlabel('Frequency, Hz','interpreter','latex','fontsize',9)
% end
% filename = 'SPL_at_NPR4p3_part2'; 
% saveas(gca, fullfile(root, filename), 'png');
%% Saving data
matfile = fullfile(root, 'OSU4p0');
save(matfile,'Fdom','SPL',"OASPL");

%%



