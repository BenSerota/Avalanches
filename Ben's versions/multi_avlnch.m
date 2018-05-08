% returns alpha-sigma plots, over multiple runs of avalanche analyses.
clear all
close all
clc
start_ben
global out_paths conds subconds num lim z_flag avl_outpath param_rows params mult_count grp cond_flag %#ok<NUSED>
DOC_basic
Avlnch_noHB_param


[Sigmas, Alphas] = deal(nan(length(tb_size),length(thresh)));

run_count = 1;
% BigStuff = cell(1,length(thresh));

if multi_flag
    for mult_count = 1:5 %length(thresh)
        
        temp = Avlnch_noHB_Gen(multi_flag,0,0); % 4 = cond 4 = ctrl
        temp = cellfun(@(x) mean(x,1), temp, 'uniformoutput', false);
        
        %% save alphas and sigmas:
        Sigmas(:,mult_count) = temp{1}; % naturally transposes
        Alphas(:,mult_count) = temp{2};
    end
else
    ...
end
%% saving
cd(avlnch_rslts)
SaveUnique('Multi_')

PlotMultiAvlnch
