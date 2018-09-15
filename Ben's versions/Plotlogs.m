function [] = Plotlogs (AVPRMS,GRP,COND,BINS)
DOC_basic

PARAM = 'alphs'; %,'ls','fls'};

if ~~COND
    condish = COND;
else
    condish = 1:4;
end

figure

for i = 1:length(condish)
    avprms = AVPRMS(GRP==condish(i));
    L = length(avprms);
    
    datasum = sparse(0);
    for ii = 1:L
        data = avprms(ii).(PARAM){5, 1};
        datasum = datasum + data;
    end
    
    %% Plotting
    len_avlnch = numel(datasum);
    probs = datasum/sum(datasum);
    series = 1:len_avlnch;
    
    %smoothing, i.e., binning:
    [bh,xs] = hist2smooth(probs,[],[],BINS);
    
    %saving for graph borders:
    BH(i) = bh(end); XS(i) = xs(end);
    
    %     loglog(xs,bh)
    plot(xs,bh,'linewidth',2)
    hold on
    % not smoothing:
    %     figure
    %     loglog(series',probs)
    
    % adding exponent line
    if i == length(condish) % end of loop
        expseries = series.^(-3/2);
        firstprob(i) = probs(1);
        plot(series,10^(log10(probs(1)))*expseries,'r--','linewidth',2) % 10^(log10(probs(1))) positions line
    end
    
    % titles
    name = 'Sizes';
    if ~~COND
        title4plot = sprintf('Probability of Neural Avalanch %s \n in %s Condition',name,conds{condish(i)});
    else
        title4plot = sprintf('Probability of Neural Avalanch %s \n in ALL Conditions',name);
    end
    
    title(title4plot,'fontsize',14,'FontWeight','bold')
    xlabel(['Avalanche' name], 'fontsize',14,'FontWeight','bold')
    ylabel('Probability','fontsize',14,'FontWeight','bold')
    
end
BH = min(BH); XS = max(XS);
set(gca, 'YScale', 'log','XScale', 'log', 'xlim',[1 XS], 'ylim',[BH 1])

% legend
condss = conds(condish);
condss{length(condish)+1} = '-3/2 exponent';
legend(condss)
end

% old:
% for i = 1:length(avprms_params)

%     %% finding last nnz element in histogrm
%     temp = flip(datasum);
%     ind = numel(datasum) - min(find((temp~=0))) + 1;
%     series = 1:ind;

%
%     switch i
%         case 1
%             name = 'Size';
%         case 2
%             name = 'Length';
%         case 3
%             name = 'Size by Length';
%     end