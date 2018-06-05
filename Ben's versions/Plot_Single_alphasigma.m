% PlotSingleAvlnch
% prepare data: MEANS and STDs + STEs
Chosen_ave_Sigmas = cellfun(@(x) mean(x),alph_sig{1},'uniformoutput',0);
Chosen_ave_Alphas = cellfun(@(x) mean(x),alph_sig{2},'uniformoutput',0);
Chosen_STD_Sigmas = cellfun(@(x) std(x),alph_sig{1},'uniformoutput',0);
Chosen_STD_Alphas = cellfun(@(x) std(x),alph_sig{2},'uniformoutput',0);
Chosen_STE_Sigmas = cell2mat(Chosen_STD_Sigmas) ./ sqrt(cell2mat(cellfun(@(x) numel(x),alph_sig{1},'uniformoutput',0)));
Chosen_STE_Alphas = cell2mat(Chosen_STD_Alphas) ./ sqrt(cell2mat(cellfun(@(x) numel(x),alph_sig{2},'uniformoutput',0)));

figure()
hold on
errorbar(cell2mat(Chosen_ave_Sigmas),cell2mat(Chosen_ave_Alphas),Chosen_STE_Sigmas,'linewidth',3)
errorbar(cell2mat(Chosen_ave_Sigmas),cell2mat(Chosen_ave_Alphas),Chosen_STE_Alphas,'horizontal','linewidth',3)
plot(cell2mat(Chosen_ave_Sigmas),cell2mat(Chosen_ave_Alphas),'-o','linewidth',2)


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

title(sprintf('Avalanche Analysis on ALL consciousness groups \n Threshold = %g STD, Time-bin = %g',thresh(chosen_th),tb_size(chosen_tb)))

xlabel('\sigma', 'fontsize',18,'FontWeight','bold')
ylabel('\alpha','fontsize',18,'FontWeight','bold')