% returns alpha-sigma plots, over multiple runs of avalanche analyses.
clear all
close all
clc
start_ben
global out_paths conds subconds thresh tb_size num lim z_flag avl_outpath param_rows params mult_count grp cond_flag pos %#ok<NUSED>
DOC_basic
Avlnch_noHB_param

% number of thresholds is either 5 (default) or more. setting here: 
if length(thresh) == 1
    l_thresh = 5;
else 
    l_thresh = length(thresh);
end

[Sigmas, Alphas] = deal(nan(l_thresh,length(tb_size)));
run_count = 1;

if multi_flag
    for mult_count = 1:l_thresh % over thresholds!
        
        temp = Avlnch_noHB_Gen(multi_flag,0,0); % 4 = cond 4 = ctrl
        temp = cellfun(@(x) mean(x,1), temp, 'uniformoutput', false);
        
        %% save alphas and sigmas:
        Sigmas(:,mult_count) = temp{1}; % naturally transposes
        Alphas(:,mult_count) = temp{2};
        fprintf('\n \n ran %g flops of avalanche analysis, out of %g \n \n', mult_count, l_thresh)
    end
else
    ...
end
%% saving
cd(avlnch_rslts)
SaveUnique('Multi_')

PlotMultiAvlnch
