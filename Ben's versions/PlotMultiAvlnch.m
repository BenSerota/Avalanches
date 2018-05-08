% PlotMultiAvlnch

f = figure('name','Multi Avalanche Analysis');
% plot(t,'-d','linewidth',3)
% hold on
for i = 1:size(Sigmas,2) 
figure(f)
scatter(Sigmas(:,i),Alphas(:,i))
plot(Sigmas(:,i),Alphas(:,i),'-d','linewidth',3)
% plot3(Sigmas(:,i),Alphas(:,i),1:10,'-d','linewidth',3) % 10 Time bins, 11 thresholds.

hold on
end

%% adding reference lines
xlim([0 3])
ylim([-3 0])
xxlim = get(gca,'xlim');
yylim = get(gca,'ylim');
plot(xxlim,[-1.5 -1.5],'black')
hold on
plot([1 1],yylim,'black')

%% aesthetics
title('Multi Avalanche Analysis')
% xlim([0 3])
% ylim([-3 0])
xlabel('\sigma', 'fontsize',18,'FontWeight','bold')
ylabel('\alpha','fontsize',18,'FontWeight','bold')
zlabel('Time Bin')

threshlegend  = cell(1,length(thresh));
for i = 1:length(thresh)
    threshlegend {i} = ['thresh = ' num2str(thresh(i))];
end
legend(threshlegend)