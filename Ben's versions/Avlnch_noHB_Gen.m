%LZC_noHB_Gen

clear all 
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
amnt_crpt = 0;
for i = 1:length(avprms)
    amnt_crpt = amnt_crpt + avprms(i).corrupt;
end
if amnt_crpt>0
    error('some rows are curropt');
else
    cd(avlnch_rslts)
    saveunique('WS_Avl')
end
%% Generate statistics
if ~exist('sts','var')  % check if doc_avprms2sts has already run
    [sts,nms,tsk,grp,subj,hs,phsps] = doc_avprms2sts(avprms,'all',1,5);
end

%% plotting
scatter_avl

%% testing significance




