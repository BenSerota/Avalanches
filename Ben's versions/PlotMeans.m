
% plot means of LC groups avalanche analyses

DOC_basic
files = {'Multi_STS_VS__2018_5_28_10_10_39',...
    'Multi_STS_MC__2018_5_28_12_8_54',...
    'Multi_STS_EMC__2018_5_28_12_55_49',...
    'Multi_STS_CTRL__2018_5_28_9_51_18'};

M_Alphas = nan(4,10);
M_Sigmas = nan(4,10);
% ~~~~~~~~~~~~~~~~~~~~~~~

for i = 1:4
    load(files{i},'Alphas'); load(files{i},'Sigmas');
    tempA = mean(Alphas,1); tempS = mean(Sigmas,1);
    
    M_Alphas(i,:) = tempA;
    M_Sigmas(i,:) = tempS;
    
end

plot(M_Sigmas',M_Alphas','-d','linewidth',3)

% adding reference lines
xlim([0 2])
ylim([-2 0])
xxlim = get(gca,'xlim');
yylim = get(gca,'ylim');
hold on
plot(xxlim,[-1.5 -1.5],'black')
hold on
plot([1 1],yylim,'black')


title(sprintf('Avalanche Analysis on all Consciousness Groups\n Averaging over 3:0.1:4.5 STDs thresholds'))
xlabel('\sigma', 'fontsize',18,'FontWeight','bold')
ylabel('\alpha','fontsize',18,'FontWeight','bold')
legend(conds)



%% Averaging over Time-Bins as well

M_M_Alphas = mean(M_Alphas,2);
M_M_Sigmas = mean(M_Sigmas,2);
A_Errs = std(M_Alphas,0,2)./sqrt(size(M_Alphas,2)); % STE for Alphas (Y axis)
S_Errs = std(M_Sigmas,0,2)./sqrt(size(M_Sigmas,2)); % STE for Sigmas (X axis)
figure()
errorbar(M_M_Sigmas,M_M_Alphas,A_Errs,'linewidth',3)
hold on
errorbar(M_M_Sigmas,M_M_Alphas,S_Errs,'horizontal','linewidth',3)
plot(M_M_Sigmas,M_M_Alphas,'-o','linewidth',3)

%Addinglabels:
text(M_M_Sigmas,M_M_Alphas,conds,'Fontsize',12,'FontWeight','bold')

xlim([0 2])
ylim([-2 0])
xxlim = get(gca,'xlim');
yylim = get(gca,'ylim');
hold on
plot(xxlim,[-1.5 -1.5],'black')
hold on
plot([1 1],yylim,'black')


title(sprintf('Avalanche Analysis on all Consciousness Groups\n Averaging over thresholds AND time-bins'))
xlabel('\sigma', 'fontsize',18,'FontWeight','bold')
ylabel('\alpha','fontsize',18,'FontWeight','bold')

