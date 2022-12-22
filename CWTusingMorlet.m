%% continous Wavelet transformation using Morlet wavelet
% source1 : https://www.mathworks.com/help/wavelet/ref/cwt.html
% source 2 : https://www.youtube.com/watch?v=jnxqHcObNK4&ab_channel=ArtemKirsanov

%% Input data
clc;
clear all;
close all;
T_aq = 2 % sec. This is the time for which sample is taken
F_aq = 204800 % Acquizition frequency of the mic
N = T_aq*F_aq % Total data collected by each mics

%fullname = fullfile('C:\Users\Kaurab Gautam\OneDrive - University of Cincinnati\Desktop\MyThesis\TestData\FarField\TR1p0\TwinMajor\AatreshFF\Major\TR_1p0\NPR_2p5', 'Pressure_FF_Pos001.dat');
root = 'C:\Users\Kaurab Gautam\OneDrive - University of Cincinnati\Desktop\MyThesis\TestData\NearField\NearField10_28\TR1p0\V1\NPR2p5' ;%'C:\Users\Kaurab Gautam\OneDrive - University of Cincinnati\Desktop\MyThesis\TestData\NearField\Twin_Rectangular_V0_NF\V1\Minor\TR_1p0\NPR_3p0';

data = load([root,'\Pressure_FF_Pos001.dat']);
data(:,1) = [];
pressure = data;
N_mics = size(pressure,1) % no of mics

%% Input for strauhaul number


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
angles = [1, 2]

%% CWT part
[wt, frq]  = cwt(pre_fluc(1,:),'amor', F_aq);
% %% Ploting part
figure;
tms = (0:1: numel(pre_fluc(1,:))-1)/F_aq;
%image(tms,frq,abs(wt),"scaled")
image("XData",tms,"YData",frq,"CData",abs(wt),"CDataMapping","direct")
axis tight
xlabel("Time (s)")
ylabel("Frequency (Hz)")








%% Scalogram
figure;
SC = wscalogram('image',wt);





