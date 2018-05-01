%LZC_noHB_Gen

% clear
close all
clc
start_ben
global out_paths conds subconds num lim z_flag avl_outpath param_rows params  %#ok<NUSED>
DOC_basic
Avlnch_noHB_param

%load aavlnch analysis?
% cd(avlnch_rslts)
% load AllSubj_AvlnchPrm_thresh2to3

%% run analysis or load data
if ~exist('avprms','var')  % check if doc_avprms2sts has already run
    cd(avlnch_rslts)
    [avprms] = Avlnch_noHB(data_frac, outname);
end


%% testing resulting avlnchs stats
% amnt_crpt = nnz([avprms.corrupt]);
% if amnt_crpt>0
%     disp('some rows are corrupt');
% %     avprms(find([avprms.corrupt] == 1)) = nan?;
% end

%% Generate statistics
if ~exist('sts','var')  % check if doc_avprms2sts has already run
%     [sts,nms,tsk,grp,subj,hs,phsps] = doc_avprms2sts(avprms,'all',1,5,5);
    [alphs,sigms,taus,gams,delts,kap,gen_kp,cos,subj,cond_subj,cond,tsk] = doc_avprms2sts_single(avprms,1,10);
    sts = zeros(33,length(alphs));
    sts0 = cat(1,sigms,alphs,taus,gams,delts,kap,gen_kp);
    sts([2,4,6,8,10,31,32],:) = sts0;
    load nec_var
    
    cd(avlnch_rslts)
    SaveUnique('WS_Avl')
end


%% registering outliers
% must be here (not after any rejection of matrices)
register_ourliers

%% reject_defects
reject_defects

%% assemble defects and reject their asses


%% plotting
% didrej = 0; % already did outlier rejection?
scatter_avl
tilefigs

%% testing significance
avlnch_var_test





