%LZC_noHB_Gen

% clear all 
close all
clc
start_ben
global out_paths conds subconds num lim z_flag avl_outpath %#ok<NUSED>
DOC_basic
Avlnch_noHB_param

%% run analysis or load data
if ~exist('avprms','var')  % check if doc_avprms2sts has already run
    cd(avlnch_rslts)
    [avprms] = Avlnch_noHB(data_frac, outname);
end
%% testing resulting avlnchs stats
amnt_crpt = nnz([avprms.corrupt]);

if amnt_crpt>0
    disp('some rows are curropt');
    avprms(find([avprms.corrupt] == 1)) = [];
end

%% Generate statistics
if ~exist('sts','var')  % check if doc_avprms2sts has already run
    [sts,nms,tsk,grp,subj,hs,phsps] = doc_avprms2sts(avprms,'all',1,5);
end

cd(avlnch_rslts)
SaveUnique('WS_Avl')

%% plotting
scatter_avl

% registering outliers
register_ourliers

%% testing significance
avlnch_ver_test

%% making new variable for 1 tb
% prmtrs = {'alphs', 'bprm1', 'bprm2', 'ls','fls','iai'};
% for i = 1:length(avprms1)
%     for ii = 1:length(prmtrs)
%         avprms1(i).(prmtrs{ii})(:,2:10) = [];
%     end
% end
[alphs,sigms,taus,gams,delts,kap,gen_kp,cos,subj,cond_subj,cond,tsk] = doc_avprms2sts_single(avprms,1,11);


