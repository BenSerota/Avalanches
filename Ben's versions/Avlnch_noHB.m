function [ALL_BIG5s] = Avlnch_noHB(data_frac,task_flag)
% this function:
% 1. takes INPUT ratio of jaco data and removes its HB component.
% 2. calculates Lempel Ziv Complexity, Zhang implementaiton.
% * task_flag = take into consideration LDGD / LDGS etc. tasks.

global out_paths subconds num lim E_T_flag %#ok<NUSED>

%% handle input
if data_frac <= 0 || data_frac > 1
    error('data ratio must be between 0 and 1')
end

if task_flag ~= 0 && task_flag ~= 1
    error('task_flag must be either 1 or 0, i.e., taking into consideration experiment task or not')
end

%% prepare
DOC_Basic2; av_param_values;
global tb_size thresh pos cont

ALL_BIG5s = cell(1,length(conds));

%% handle data_frac < 1
if data_frac < 1
    small_amnt_sbjcts = round(data_frac * amnt_sbjcts); %#ok
    % handle case where # subj is 0 due to low data_frac
    if nnz(small_amnt_sbjcts==0)>0
        error('data_frac appears to be too low, number of subjects in some conditions is 0')
    end
    amnt_sbjcts = num2cell(amnt_sbjcts);
    small_amnt_sbjcts = num2cell(small_amnt_sbjcts);    % moving to cell structure for @cellfun only
    names_inds = cellfun(@(x,y) randperm(x,y), amnt_sbjcts, small_amnt_sbjcts,'UniformOutput' ,false);
    
    % choosing subjects
    temp = cell(length(names_inds),1);
    for i = 1:length(names_inds)
        temp{i} = NAMES{i}(names_inds{i}); %#ok
    end
    
    NAMES = temp; clear temp;
    small_amnt_sbjcts = cell2mat(small_amnt_sbjcts); % returning to mat structure
    amnt_sbjcts = small_amnt_sbjcts; % giving back to amnt_sbjcts for JacoClock !
end

%% go
% main outpus:
% [f,bprm,bprm2,ls,fls,iai] = deal(cell(1,length(conds)));

% % for rate = rates
%     while ~finito
%         [f{cnd}{subj}, bprm{cnd}{subj}, bprm2{cnd}{subj},...
%             ls{cnd}{subj}, fls{cnd}{subj}, iai{cnd}{subj}]...
%             = Do1Sbj(NAMES, cnd, subj); % allocate LZC scores (1 grade or 4... egal)
%         [cnd, subj, finito] = JacoClock(amnt_sbjcts, cnd, subj);    % advaning us in Jaco clock
%     end
% % end

% main outpus: BIG5s

% for rate = rates
    while ~finito
        [ALL_BIG5s{cnd}{subj}] = Do1Sbj(NAMES, cnd, subj); % allocate LZC scores (1 grade or 4... egal)
        [cnd, subj, finito] = JacoClock(amnt_sbjcts, cnd, subj);    % advaning us in Jaco clock
    end
% end

%% assisting functions


function [BIG5] = Do1Sbj(NAMES, cnd, subj) % NOTE: consider making task_flag an input
global out_paths subconds task_flag tb_size thresh pos cont

cd(out_paths{cnd})
name = char(NAMES{cnd}(subj));
name_p = [ name '_prep'];
name_I = [ name '_HBICs'];

% load each set of ICs per subj
load(name_p);                                                              % load subject
load(name_I);                                                              % load HB info

%% TO DO :
%TODO: here a big mss to fix/write: resave in WS data with no hb. check
% then matcat it 
% then build BIG5?
%%

for i  = 1:length(subconds)
     DATAnhb{i} = RemoveHB1Task(i,final,Comps2Reject);
     DATAnhb{i} = DATAnhb{i}'; % cuz tomer's code wants mat to be : datapoints X channels.
     [BIG5] = eeg2avalnch_ben(DATAnhb{i},tb_size,thresh,pos,cont)
end


%%
% temp = data;
% temp_raw = matcat(temp);
% raw_data{ii} = temp_raw'; % must have correct structure!
% clear temp
% BIG5 = eeg2avalnch_ben(raw_data{ii},tb_size,thresh,pos,cont);
% SUBJECTS{ii} = BIG5; % feeding into main cell array

clear final Comps2Reject   % just in case...

function [DATAnhb] = RemoveHB1Task(i,final,Comps2Reject)
global subconds num lim
DATAwhb.data = final.(subconds{i}).data;
DATAwhb.w = final.(subconds{i}).w;
DATAwhb.sph = final.(subconds{i}).sph;
rejcomps = Comps2Reject.(subconds{i});

DATAnhb = remove_HB(DATAwhb, rejcomps, num, lim);

clear DATAwhb rejcomps

DATAnhb = reshape(DATAnhb,206,385,[]);



function [DATAnhb] = remove_HB(DATAwhb, rejcomps, num, lim)
% global num lim 
% This code is meant to remove a number of hb ICs form our data.

% ICA activation matrix  = W * Sph * Raw:
ICAact = DATAwhb.w * DATAwhb.sph * DATAwhb.data;

% choosing #num comps from all hb cmps
if isempty(rejcomps)
    DATAnhb = DATAwhb.data; % no need to eliminate any comps
    return
elseif num > length(rejcomps)   % in case we are asking for too many comps than there are left
    num = length(rejcomps);     % num = number of hb comps
end

% treating only comps below limit
rejcomps(rejcomps>lim) = [];
if isempty(rejcomps)    % if we are left with no comps, stop
    DATAnhb = DATAwhb.data;
    return
elseif num > length(rejcomps) % in case we are asking for too many comps than there are left (after >lim elimination)
    num = length(rejcomps);
end

rejcomps = rejcomps(1:num); % take first 'num comps

%elimninating hb
ICAact(rejcomps,:) = [];
w_inv = pinv(DATAwhb.w*DATAwhb.sph);
w_inv(:,rejcomps) = [];  % rejecting corresponding colomns
DATAnhb = w_inv * ICAact;


function [cnd,subj, finito] = JacoClock(amnt_sbjcts, cnd, subj)

subj = subj + 1;
if subj > amnt_sbjcts(cnd)
    cnd = cnd + 1;
    subj = 1;
end

finito = 0;

if cnd > 4
    finito = 1;
end

% function [] = SaveUniqueName(root_name)
% if ~isstring(root_name) && ~ischar(root_name)
%     error('input must be of class char or string')
% end
% cl = fix(clock);
% stamp = strrep(mat2str(cl),' ','_');
% stamp = strrep(stamp,'[','');
% stamp = strrep(stamp,']','');
% UniqueName = [root_name '_' stamp];
% cd ('E:\Dropbox\Ben Serota\momentary\WS')
% save (UniqueName)
