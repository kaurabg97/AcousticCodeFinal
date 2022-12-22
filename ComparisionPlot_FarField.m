clear all;
close all;

root = 'C:\Users\Kaurab Gautam\OneDrive - University of Cincinnati\Desktop\MyThesis';
dataA= load([root,'\ReducedShroud3p6.mat'])
dataB= load([root,'\NoShroud3p6.mat'])
%dataC= load([root,'\OSUAA90A1p55.mat'])
Fdom= load([root,'\Fdom.mat'])
angles = [152, 148, 144, 140, 136, 132, 126, 120, 116, 110, 105,100, 90, 70, 60, 45]
%% Plotting
txt = input("Enter a value: " )
if txt == 1 % 1 fo5 without cover2 2 for cover on flange

    n_mics =16
    figure;
    plot(angles,dataA.OASPL,"R", 'LineWidth',1.1);
    hold on;
    plot(angles,dataB.OASPL,"b", 'LineWidth',1.1);
    hold on;
    % plot(angles,dataC.OASPL_rms,"g", 'LineWidth',1.1)
    % hold on;
   % plot(angles,dataC.OASPL_rms,"k", 'LineWidth',1.1);
   % hold on;
    % plot(angles,dataD.OASPL_rms,"m", 'LineWidth',1.1)
    % hold on
    xlim([45 155]);
    ylim([105 135]);
    grid on;
    legend('Reduced shroud','No shroud','Shroud','interpreter','latex','fontsize',11);
    xlim([45 155]);
    ylim([105 140]);
    txt = [' OASPL plot for NPR 3.6'];
    title(txt, 'interpreter','latex','fontsize',11);
    ylabel('OASPL, dB','interpreter','latex','fontsize',11);
    xlabel(' Observation angles ( $ \psi $ ), Hz','interpreter','latex','fontsize',11);

    %% SPL Plots
    figure;
    subplot(1,3,1)
    semilogx(Fdom.Fdom,dataA.SPL(:,1),"R")
    hold on
    semilogx(Fdom.Fdom,dataB.SPL(:,1),'b')
    hold on
    % semilogx(Fdom.Fdom,dataC.SPL(:,1),'g')
    % hold on
    %semilogx(Fdom.Fdom,dataC.SPL(:,1),'k')
    % hold on
    % semilogx(Fdom.Fdom,dataE.SPL(:,1),'m')
    xlim([300 150000])
    ylim([60 125])
    grid on
    txt = [' $ \psi $ =', num2str(152), ' degree']
    title(txt, 'interpreter','latex','fontsize',11)
    ylabel('SPL, dB','interpreter','latex','fontsize',11)
    xlabel('Frequency, Hz','interpreter','latex','fontsize',11)
    legend('Reduced shroud','No Shroud','interpreter','latex','fontsize',11)
    %legend('Aatresh Test','Baseline','Baseline (No Cover)','Shroud','Shroud ( No cover')
    pbaspect([1 1 1])

    subplot(1,3,2)
    semilogx(Fdom.Fdom,dataA.SPL(:,13)+10,"R")
    hold on
    semilogx(Fdom.Fdom,dataB.SPL(:,13)+20,'b')
    hold on
    % semilogx(Fdom.Fdom,dataC.SPL(:,13),'g')
    % hold on
%     semilogx(Fdom.Fdom,dataC.SPL(:,13),'k')
%     hold on
    % semilogx(Fdom.Fdom,dataE.SPL(:,13),'m')
    xlim([300 150000])
    ylim([60 125])
    grid on
    txt = [' $ \psi $ =', num2str(90), ' degree']
    title(txt, 'interpreter','latex','fontsize',11)
    ylabel('SPL, dB','interpreter','latex','fontsize',11)
    xlabel('Frequency, Hz','interpreter','latex','fontsize',11)
    legend('Reduced shroud','No Shroud','interpreter','latex','fontsize',11)
    pbaspect([1 1 1])

    subplot(1,3,3)
    semilogx(Fdom.Fdom,dataA.SPL(:,16)+10,"R")
    hold on
    semilogx(Fdom.Fdom,dataB.SPL(:,16)+20,'b')
    hold on
    % semilogx(Fdom.Fdom,dataC.SPL(:,16),'g')
    % hold on
%     semilogx(Fdom.Fdom,dataC.SPL(:,16),'k')
%     hold on
    % semilogx(Fdom.Fdom,dataE.SPL(:,16),'m')
    xlim([300 150000])
    ylim([60 125])
    grid on
    txt = [' $ \psi $ =', num2str(45), ' degree']
    title(txt, 'interpreter','latex','fontsize',11)
    ylabel('SPL, dB','interpreter','latex','fontsize',11)
    xlabel('Frequency, Hz','interpreter','latex','fontsize',11)
    legend('Reduced shroud','No Shroud','interpreter','latex','fontsize',11)
    pbaspect([1 1 1])
else
    n_mics =16
    figure;
    plot(angles,dataA.OASPL_rms,"R", 'LineWidth',1.1)
    hold on;
    %     plot(angles,dataB.OASPL_rms,"b", 'LineWidth',1.1)
    %     hold on;
    plot(angles,dataC.OASPL_rms,"b", 'LineWidth',1.1)
    hold on;
    %     plot(angles,dataD.OASPL_rms,"k", 'LineWidth',1.1)
    %     hold on;
    plot(angles,dataC.OASPL_rms,"k", 'LineWidth',1.1)
    hold on
    xlim([45 155])
    ylim([105 135])
    grid on
    legend('Aatresh Test','Baseline','Shroud','interpreter','latex','fontsize',11)
    xlim([45 155])
    ylim([105 140])
    txt = [' OASPL plot for NPR 3.67 (Flange uncovered)']
    title(txt, 'interpreter','latex','fontsize',11)
    ylabel('OASPL, dB','interpreter','latex','fontsize',11)
    xlabel(' Observation angles ( $ \psi $ ), Hz','interpreter','latex','fontsize',11)

    %% SPL Plots
    figure;
    subplot(1,3,1)
    semilogx(Fdom.Fdom,dataA.SPL(:,1),"R")
    hold on
%     semilogx(Fdom.Fdom,dataB.SPL(:,1),'b')
%     hold on
    semilogx(Fdom.Fdom,dataC.SPL(:,1),'g')
    hold on
    %semilogx(Fdom.Fdom,dataD.SPL(:,1),'k')
    % hold on
     semilogx(Fdom.Fdom,dataE.SPL(:,1),'m')
    xlim([300 150000])
    ylim([60 125])
    grid on
    txt = [' $ \psi $ =', num2str(152), ' degree']
    title(txt, 'interpreter','latex','fontsize',11)
    ylabel('SPL, dB','interpreter','latex','fontsize',11)
    xlabel('Frequency, Hz','interpreter','latex','fontsize',11)
    legend('Aatresh Test','Baseline','Shroud','interpreter','latex','fontsize',11)
    %legend('Aatresh Test','Baseline','Baseline (No Cover)','Shroud','Shroud ( No cover')
    pbaspect([1 1 1])

    subplot(1,3,2)
    semilogx(Fdom.Fdom,dataA.SPL(:,13),"R")
    hold on
    semilogx(Fdom.Fdom,dataB.SPL(:,13),'b')
    hold on
    % semilogx(Fdom.Fdom,dataC.SPL(:,13),'g')
    % hold on
    semilogx(Fdom.Fdom,dataC.SPL(:,13),'k')
    hold on
    % semilogx(Fdom.Fdom,dataE.SPL(:,13),'m')
    xlim([300 150000])
    ylim([60 125])
    grid on
    txt = [' $ \psi $ =', num2str(90), ' degree']
    title(txt, 'interpreter','latex','fontsize',11)
    ylabel('SPL, dB','interpreter','latex','fontsize',11)
    xlabel('Frequency, Hz','interpreter','latex','fontsize',11)
    legend('Aatresh Test','Baseline','Shroud','interpreter','latex','fontsize',11)
    pbaspect([1 1 1])

    subplot(1,3,3)
    semilogx(Fdom.Fdom,dataA.SPL(:,16),"R")
    hold on
    semilogx(Fdom.Fdom,dataB.SPL(:,16),'b')
    hold on
    % semilogx(Fdom.Fdom,dataC.SPL(:,16),'g')
    % hold on
    semilogx(Fdom.Fdom,dataC.SPL(:,16),'k')
    hold on
    % semilogx(Fdom.Fdom,dataE.SPL(:,16),'m')
    xlim([300 150000])
    ylim([60 125])
    grid on
    txt = [' $ \psi $ =', num2str(45), ' degree']
    title(txt, 'interpreter','latex','fontsize',11)
    ylabel('SPL, dB','interpreter','latex','fontsize',11)
    xlabel('Frequency, Hz','interpreter','latex','fontsize',11)
    legend('Aatresh Test','Baseline','Shroud')
    pbaspect([1 1 1])
end




% for i = 1:n_mics/2
%     ax(i)= subplot(2,4,i)
%     semilogx(Fdom.Fdom,dataA.SPL(:,i),"R")
%     hold on
%     semilogx(Fdom.Fdom,dataB.SPL(:,i),'b')
%     hold on
%     semilogx(Fdom.Fdom,dataC.SPL(:,i),'g')
%     hold on
%     semilogx(Fdom.Fdom,dataD.SPL(:,i),'k')
%     hold on
%     semilogx(Fdom.Fdom,dataE.SPL(:,i),'m')
%     hold on
%     xlim([300 150000])
%     ylim([60 125])
%     grid on
%     txt = [' $ \psi $ =', num2str(angles(i)), ' degree']
%     title(txt, 'interpreter','latex','fontsize',11)
%     ylabel('SPL, dB','interpreter','latex','fontsize',11)
%     xlabel('Frequency, Hz','interpreter','latex','fontsize',11)
%     legend('Aatresh Test','Baseline','Baseline (No Cover)','Shroud','Shroud ( No cover')
% end
% figure;
%
% for i = 1:n_mics/2
%     ax(i)= subplot(2,4,i)
%     ax(i)= subplot(2,4,i)
%     semilogx(Fdom.Fdom,dataA.SPL(:,i+(n_mics/2)),"R")
%     hold on
%     semilogx(Fdom.Fdom,dataB.SPL(:,i+(n_mics/2)),'b')
%     hold on
%     semilogx(Fdom.Fdom,dataC.SPL(:,i+(n_mics/2)),'g')
%     hold on
%     semilogx(Fdom.Fdom,dataD.SPL(:,i+(n_mics/2)),'k')
%     hold on
%     semilogx(Fdom.Fdom,dataE.SPL(:,i+(n_mics/2)),'m')
%     hold on
%     grid on
%     txt = [' $ \psi $ =', num2str(angles(i+(n_mics/2))), ' degree']
%     title(txt, 'interpreter','latex','fontsize',11)
%     ylabel('SPL, dB','interpreter','latex','fontsize',11)
%     xlabel('Frequency, Hz','interpreter','latex','fontsize',11)
%     legend('Aatresh Test','Baseline','Baseline (No Cover)','Shroud','Shroud ( No cover')
% end

