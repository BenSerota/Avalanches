% PlotMultiAvlnch

f = figure('name','Multi Avalanche Analysis');
g = figure('name','Multi Avalanche Analysis');

for ii = [f g]
    figure(ii)
    if ii == f
        plot(Sigmas',Alphas','-d','linewidth',3)
    else
        plot(mean(Sigmas,1),mean(Alphas,1),'-d','linewidth',3)
    end
    %% adding reference lines
    xlim([0 2])
    ylim([-2 0])
    xxlim = get(gca,'xlim');
    yylim = get(gca,'ylim');
    hold on
    plot(xxlim,[-1.5 -1.5],'black')
    hold on
    plot([1 1],yylim,'black')
    
    %% aesthetics
    if l_thresh ~= 1
        thresh4title = [thresh(1), thresh(2)-thresh(1), thresh(numel(success))];
        thresh4title = [num2str(thresh4title(1)) ':' num2str(thresh4title(2)) ':' num2str(thresh4title(3))];
    else
        thresh4title = num2str(round(thresh,1));
        thresh4title = [thresh4title ' +- 0.2 Prcntiles'];
    end
    
    if ii == f
        title(sprintf('Avalanche Analysis on %s group \n thresh = %s STDs',out_b,thresh4title))
    else
        title(sprintf('Avalanche Analysis on %s group \n AVERAGE over thresholds: %s STDs' ,out_b,thresh4title))
    end
    
    xlabel('\sigma', 'fontsize',18,'FontWeight','bold')
    ylabel('\alpha','fontsize',18,'FontWeight','bold')
    
    threshlegend  = cell(1,length(success)); % using success and not l_thresh
    if length(thresh) == 1 % in case thresh is stated in single term (not vector)
        for i = 1:length(success)
            s = thresh + (-3 + i)*0.1;
            threshlegend {i} = sprintf('top %g percentile', s);
        end
    else
        if ii == f % no need for legend when averaging
            for i = 1:length(success)
                threshlegend {i} = sprintf('thresh = %g STDs', thresh(i));
            end
            legend(threshlegend)
        end
    end
    
    
end