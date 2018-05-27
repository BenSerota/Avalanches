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
    if length(thresh) ~= 1
        thresh = round(thresh);
        thresh4title = [thresh(1), thresh(2)-thresh(1), thresh(end)];
        thresh4title = num2str(thresh4title);
    else
        thresh4title = num2str(thresh);
        thresh4title = [thresh4title ' +- 0.1 STDs'];
    end
    
    if ii == f
        title(sprintf('Avalanche Analysis on %s group, thresh = %s, NOT averaging over thresholds',out_b,thresh4title))
    else
        title(sprintf('Avalanche Analysis on %s group, thresh = %s, AVERAGE over thresholds' ,out_b,thresh4title))
    end
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
    
end