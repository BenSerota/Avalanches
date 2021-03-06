function [sig_alph,sts] = LZC_noHB_Gen (multi_flag, var_flag, scat_flag)
% runs avalanch analysis over JAco DOC database
% multi_flag = are we creating a alpha-sigma plot?
% cond_flag (globally inserted) = analyze statistics for condition # or all(0) ?
% var_flag = should we also perform ANOVA?
% scat_flag = should we scatter individual grades? (per dataset per parammeter).

close all
start_ben
global out_paths conds subconds tb_size num lim z_flag avlnch_rslts param_rows params mult_count grp cond_flag out_b %#ok<NUSED>
DOC_basic
Avlnch_noHB_param

%% load avalnche analysis?
cd(avlnch_rslts)
%      load AllSubj_AvlnchPrm_thresh2to3
%      load AllSubj_AvlnchPrm_thresh3to4
%      load AllSubj_AvlnchPrm_thresh005
%      load AllSubj_AvlnchPrm_thresh001
%      load('NotAve_Cond = CTRL_thresh_3_5')
load NotAve_CTRL_thresh_3_4p5
load group_list

%% run analysis or load data
if ~exist('avprms','var')  % check if doc_avprms2sts has already run
    cd(avlnch_rslts)
    [avprms] = Avlnch_noHB(data_frac, outname);
end


%% testing resulting avlnchs stats
if ~multi_flag || (multi_flag && mult_count == 1) % if running multiple, do this only once
    amnt_crpt = nnz([avprms.corrupt]);
    if amnt_crpt>0
        fprintf('\n %g datasets are corrupt \n',amnt_crpt);
    else
        fprintf('\n GOOD JOB there Dude! no datasets are corrupt \n');
    end
    
    crpt = find([avprms.corrupt] == 1)';
    
    avprms = avprms([avprms.corrupt] == 0,:); % throwing away bad rows in avprms
    grp(crpt) = []; % throwing away bad rows in grp
    amnt_sbjcts_analyzed = [...
        length(unique([avprms(find([avprms.cond] == 1)).cond_subj])),...
        length(unique([avprms(find([avprms.cond] == 2)).cond_subj])),...
        length(unique([avprms(find([avprms.cond] == 3)).cond_subj])),...
        length(unique([avprms(find([avprms.cond] == 4)).cond_subj])) ];
    amnt_dsets_analyzed = [nnz(grp == 1),nnz(grp == 2),nnz(grp == 3),nnz(grp == 4)];
end

%% Generate statistics
if ~exist('sts','var')  % check if doc_avprms2sts has already run
    if cond_flag == 1 || cond_flag == 2 || cond_flag == 3 || cond_flag == 4
        avprms = avprms(grp == cond_flag);
    elseif cond_flag ~= 0
        error('cond_flag must be an integer between 0 and 4)')
    end
    
    if multi_flag
        [sts,nms,tsk,grp,subj,hs,phsps] = doc_avprms2sts(avprms,'all',tb_size(1),tb_size(end),1,0,1,mult_count,mult_count);
    else
        % for single run
        [alphs,sigms,taus,gams,delts,kap,gen_kp,cos,subj,cond_subj,cond,tsk] = doc_avprms2sts_single(avprms,chosen_tb,chosen_th);
        sts = zeros(33,length(alphs));
        sts0 = cat(1,sigms,alphs,taus,gams,delts,kap,gen_kp,cos);
        sts(param_rows,:) = sts0;
        load nec_var2
    end
end


%% registering outliers
% must be here (not after any rejection of matrices)
register_ourliers

%% reject_defects
% assemble defects and reject their asses
reject_defects

if multi_flag
    for i = 1:length(phsps)
        phsps{i}(bads,:) = [];
    end
    sig_alph{1} = phsps{1};
    sig_alph{2} = phsps{2};
    sts = [];
else
    sig_alph = [];
end

%% plotting
if scat_flag
    scatter_avl
    tilefigs
end
%% testing significance
if var_flag
%     avlnch_var_test
    Avl_vartest
end

Plotlogs(avprms,grp,0,bins)
Plotlogs(avprms,grp,1,bins)
Plotlogs(avprms,grp,2,bins)
Plotlogs(avprms,grp,3,bins)
Plotlogs(avprms,grp,4,bins)
tilefigs



