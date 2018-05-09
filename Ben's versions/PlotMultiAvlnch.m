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
xlim([0 2])
ylim([-2 0])
xxlim = get(gca,'xlim');
yylim = get(gca,'ylim');
plot(xxlim,[-1.5 -1.5],'black')
hold on
plot([1 1],yylim,'black')

%% aesthetics
if length(thresh) ~= 1
    thresh = round(thresh);
    thresh4title = [thresh(1), thresh(2)-thresh(1), thresh(end)];
    thresh4title = num2str(thresh4title);
else
    thresh4title = num2str(thresh);
    thresh4title = [thresh4title ' +- 0.1 STDs'];
end

title(['Multi Avalanche Analysis, thresh = ' thresh4title])
% xlim([0 3])
% ylim([-3 0])
xlabel('\sigma', 'fontsize',18,'FontWeight','bold')
ylabel('\alpha','fontsize',18,'FontWeight','bold')
zlabel('Time Bin')

threshlegend  = cell(1,l_thresh);
if length(thresh) == 1 % in case thresh is stated in single term (not vector)
    for i = 1:l_thresh
        s = (-3 + i)*0.1;
        if s<0
            s = [ num2str(s) 'STD'];
        elseif s>0
            s = ['+' num2str(s) 'STD'];
        elseif s ==0
            s = [];
        end
        
        threshlegend {i} = sprintf('top %g percentile %s', thresh,s);
        
    end
else
    
    
    for i = 1:l_thresh
        threshlegend {i} = sprintf('thresh = %g', thresh(i));
    end
end
legend(threshlegend)