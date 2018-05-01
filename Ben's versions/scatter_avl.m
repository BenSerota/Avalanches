
%% plotting

%% plot
scatter_prms_by_cond(sts,grp,ext);
cd(avlnch_rslts)


%% inner functions
function [] = scatter_prms_by_cond(prms,grp,ext)
global conds param_rows params

for i = 1:length(params) % over parameters
    figure('name',params{i})
    row = param_rows(i);  % parameter row
    for ii = 1:length(conds) %conds
        toscatter = prms(row , grp == ii);
        scatter(ii*ones(1,length(toscatter)),toscatter)
        means(ii) = mean(toscatter);
        hold on
        scatter(ii,means(ii),'bd')
        extremes = toscatter(ext(row,grp==ii));
        scatter(ii*ones(1,length(extremes)),extremes,'k*')
        
        % setting figure borders
        uplim(ii) = max(toscatter);
        lowlim(ii) = min(toscatter);
    end
    plot(1:4,means)
    title([params{i} ' per lvl of Consciousness'])
    xticks(1:4)
    set(gca,'xticklabel',conds)
    xlabel('Level of Consciousness')
    xlim([0.5 4.5])
    uplim = max(uplim);
    lowlim = min(lowlim);
    ylim([lowlim-.2 uplim+.2]) %
end
tilefigs
end

