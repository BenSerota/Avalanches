function [] = plotlier(cnd, subj_in_cond, task, dispchan, winlen)
%  plotlier(ondition, subj_in_cond, task, # channels to display, time window size)

% from LZC_noHB_Gen
% clear
close all
start_ben
global out_paths conds subconds num lim z_flag avl_outpath param_rows params %#ok<NUSED>
DOC_basic
Avlnch_noHB_param

%load aavlnch analysis?
% cd(avlnch_rslts)
% load AllSubj_AvlnchPrm_thresh2to3

%% run analysis or load data
if ~exist('avprms','var')  % check if doc_avprms2sts has already run
    cd(avlnch_rslts)
    [avprms] = Avlnch_noHB2(data_frac, outname, cnd, subj_in_cond, task, dispchan, winlen);
end

function avprms = Avlnch_noHB2(data_frac,outname, cnd, subj_in_cond, task, dispchan, winlen)
% this function:
% 1. takes INPUT ratio of jaco data and removes its HB component.
% 2. calculates Avalanches parameters.
% 3. saves output as it goes. able to continue from point at which it stopped.

global out_paths subconds num lim tb_size thresh pos z_flag avl_outpath %#ok<NUSED>

%% handle input
if data_frac <= 0 || data_frac > 1
    error('data ratio must be between 0 and 1')
end

%% prepare
DOC_basic2; %av_param_values;

%% handle data_frac < 1
if data_frac < 1
    small_amnt_sbjcts = round(data_frac * amnt_sbjcts); %#ok
    % handle case where # subj is 0 due to low data_frac
    if nnz(small_amnt_sbjcts==0)>0
        error('data_frac appears to be too low, number of subjects in some conditions is 0')
    end
    amnt_sbjcts = num2cell(amnt_sbjcts);
    small_amnt_sbjcts = num2cell(small_amnt_sbjcts);                        % moving to cell structure for @cellfun only
    names_inds = cellfun(@(x,y) randperm(x,y), amnt_sbjcts, small_amnt_sbjcts,'UniformOutput' ,false);
    
    % choosing subjects
    temp = cell(length(names_inds),1);
    for i = 1:length(names_inds)
        temp{i} = NAMES{i}(names_inds{i}); %#ok
    end
    
    NAMES = temp; clear temp;
    small_amnt_sbjcts = cell2mat(small_amnt_sbjcts);                        % returning to mat structure
    amnt_sbjcts = small_amnt_sbjcts;                                        % giving back to amnt_sbjcts for JacoClock !
end


flds = {'subj' 'cond_subj' 'cond' 'tsk' 'alphs' 'bprm1' 'bprm2' 'ls' 'fls' 'iai' 'corrupt'};
avprms = cell2struct(cell(183*4,length(flds)),flds,2);                  %183x4 = subjects x subconditions
counter = 0;

%% go
[avprms,~] = Do1Sbj2(NAMES, cnd, subj_in_cond,avprms,counter, z_flag , task, dispchan, winlen);
cd(avl_outpath)
save(outname)

%% assisting functions

function [avprms,counter] = Do1Sbj2(NAMES, cnd, subj_in_cond, avprms,counter, z_flag, task, dispchan, winlen)
global out_paths subconds tb_size thresh pos

cd(out_paths{cnd})
name = char(NAMES{cnd}(subj_in_cond));
name_p = [ name '_prep'];
name_I = [ name '_HBICs'];

% load each set of ICs per subj
load(name_p);                                                              % load subject
load(name_I);                                                              % load HB info
tic
i  = task;
teeg = RemoveHB1Task2(i,final,Comps2Reject)';
[smpls,sensr] = size(teeg);
epoch_N = smpls/385;
teeg = reshape(teeg,[385 epoch_N sensr]);
if z_flag
    teeg = zscore(teeg);
end
teeg = cat(1,teeg,zeros([1,epoch_N,sensr]));
teeg = reshape(teeg,[epoch_N*386,sensr]);

teeg = teeg';

eegplot_CPL(teeg,'srate',250,'dispchans',dispchan,'winlength',winlen);

function [DATAnhb] = RemoveHB1Task2(i,final,Comps2Reject)
global subconds num lim
DATAwhb.data = final.(subconds{i}).data;
DATAwhb.w = final.(subconds{i}).w;
DATAwhb.sph = final.(subconds{i}).sph;
rejcomps = Comps2Reject.(subconds{i});

DATAnhb = remove_HB2(DATAwhb, rejcomps, num, lim);

clear DATAwhb rejcomps % just in case...


function [DATAnhb] = remove_HB2(DATAwhb, rejcomps, num, lim)
% global num lim
% This code is meant to remove a number of hb ICs form our data.
% ICA activation matrix  = W * Sph * Raw:
ICAact = DATAwhb.w * DATAwhb.sph * DATAwhb.data;

% choosing #num comps from all hb cmps
if isempty(rejcomps)
    DATAnhb = DATAwhb.data;                                                 % i.e. no need to eliminate any comps
    return
elseif num > length(rejcomps)                                               % in case we are asking for too many comps than there are left
    num = length(rejcomps);                                                 % num = number of hb comps
end

% treating only comps below limit
rejcomps(rejcomps>lim) = [];
if isempty(rejcomps)                                                        % if we are left with no comps, stop
    DATAnhb = DATAwhb.data;
    return
elseif num > length(rejcomps)                                               % in case we are asking for too many comps than there are left (after >lim elimination)
    num = length(rejcomps);
end

rejcomps = rejcomps(1:num);                                                 % take first 'num comps

%elimninating hb
ICAact(rejcomps,:) = [];
w_inv = pinv(DATAwhb.w*DATAwhb.sph);
w_inv(:,rejcomps) = [];                                                     % rejecting corresponding colomns
DATAnhb = w_inv * ICAact;

