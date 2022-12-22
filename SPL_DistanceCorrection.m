SPL_data = load('C:\Users\Kaurab Gautam\OneDrive - University of Cincinnati\Desktop\MyThesis\TestData\OSU\050272022_FF_AA90\OSUAA90A1p30.mat')
M = [2.72, 2.82, 3.01, 3.32, 3.58, 3.78, 3.35, 3.16]
T = 1.135634; 
n_mics = size(SPL_data.SPL,2)
Scaled_SPL(size(SPL_data.SPL))= 0;
for i = 1: n_mics
    Scaled_SPL(:,i) = SPL_data.SPL(:,i)-20*log10(T/M(i));
end
Fdom = load('C:\Users\Kaurab Gautam\OneDrive - University of Cincinnati\Desktop\MyThesis\TestData\OSU\050272022_FF_AA90\Fdom.mat')
%% Plotting SPL
figure(1);
angles = [150 145, 140, 135, 130, 120, 105, 90]
for i = 1:n_mics
    %ax(i)= subplot(2,3,i)
    semilogx(Fdom.Fdom(1,:),Scaled_SPL(:,i), LineWidth=2)
    xlim([300 104800])
    %ylim([60 125])
    grid on
    txt = [' $ \psi $ =', num2str(angles(i)), ' degree'];
    title(txt, 'interpreter','latex','fontsize',9)
    ylabel('SPL, [dB]','interpreter','latex','fontsize',9)
    xlabel('St','interpreter','latex','fontsize',9)
    hold on
    legend
end
%% Savingd data
matfile = fullfile(root, 'OSU3p0');
save(matfile,'SPL',"OASPL");