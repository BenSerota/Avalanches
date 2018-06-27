
function [RHO, PVAL] = Corr_LZC_AvPrm(PRM_FLAG)
% calculates Pearson correlation between pre-calculated LZC scores and
% Alphas or Sigmas values.
% PRM_FLAG is 0/1 = sigmas / alphas.

if PRM_FLAG
    PRM = 'alphs';
else
    PRM = 'sigms';
end

load('allLZCs.mat')
load('th=3p4_tb=1','bads')
PRM_STR = load('th=3p4_tb=1',PRM);
PARAM = PRM_STR.(PRM);
clear PRM_STR
load('group_list.mat')

LZC_flats = cellfun(@(x) x', LZC, 'uniformoutput',0);
LZC_flats = cellfun(@(x) x(:), LZC_flats, 'uniformoutput',0);
LZC_flat = cat(1,LZC_flats{1},LZC_flats{2},LZC_flats{3},LZC_flats{4});
LZC_flat(bads) = [];
grp(bads) = [];
PARAM = PARAM';

% normalize
% LZC_flat = zscore(LZC_flat);
% alphs = zscore(alphs);

info = histc(grp, unique(grp)); % how much data from each cond.
n = min(info);

rep = 1e5;
% in order not to too much weight to the conds with more data:
% let's randperm each time the minimal amount of non-rej subj (48),
% correlate, and repeat 1000 times and and then average over all correlations.
% notice, choosing indeces OF GRP, for ease.

% preallocate
[rho,pval] = deal(nan(1,rep));

for i = 1:rep
    chosen_vs = randperm(info(1),n);
    chosen_mcs = randi([info(1)+1,sum(info(1:2))],[1 n]);
    chosen_emcs = randi([sum(info(1:2))+1,sum(info(1:3))],[1 n]);
    chosen_ctr = 1:info(4); % all ctr
    
    norm_lzcs = [LZC_flat(chosen_vs), LZC_flat(chosen_mcs), LZC_flat(chosen_emcs), LZC_flat(chosen_ctr)];
    norm_lzcs = norm_lzcs (:);
    norm_prm = [PARAM(chosen_vs);PARAM(chosen_mcs);PARAM(chosen_emcs);PARAM(chosen_ctr)];
    
    % correlation
    [rho(i),pval(i)] = corr(norm_lzcs,norm_prm);
    
    % for last run : plot for niceliness
    if i == rep
        p = polyfit(norm_lzcs,norm_prm,1);
        yfit = polyval(p,norm_lzcs);
        plot(norm_lzcs,yfit,'g','linewidth',2)
        
        hold on
        scatter(norm_lzcs,norm_prm)
        PRM(end)='a'; PRM(1) = upper(PRM(1)); % just for titles
        title(sprintf('Pearson Correlation \n Lempel Ziv Complexity Grade X %s value',PRM),'fontsize',14,'FontWeight','bold')
        xlabel('LZC score', 'fontsize',14,'FontWeight','bold')
        PRM(1) = lower(PRM(1));
        ylabel(sprintf('\\%s',PRM),'fontsize',20,'FontWeight','bold')
        xlim([min(norm_lzcs) 1])
    end
    
end

FISHER = mean(atanh(rho));
R = tanh(FISHER);
text(0.9,mean(norm_prm),sprintf('R = %g', round(R,2)),'fontsize',12,'FontWeight','bold')

