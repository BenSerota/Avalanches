function avprms = Avlnch_noHB(data_frac,outname)
% this function:
% 1. takes INPUT ratio of jaco data and removes its HB component.
% 2. calculates Lempel Ziv Complexity, Zhang implementaiton.
% * task_flag = take into consideration LDGD / LDGS etc. tasks.

global out_paths subconds num lim E_T_flag %#ok<NUSED>

%% handle input
if data_frac <= 0 || data_frac > 1
    error('data ratio must be between 0 and 1')
end


%% prepare
DOC_Basic2; av_param_values;


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

%
if exist('outname','var')
    try
        load(outname)
        fnams = cellstr(char(avprms.file));%#ok
        rthere = cellfun(@isempty,fnams);
        fnams(rthere) = '';
        avprms(rthere) = '';
        counter = length(avprms);
        k = length(union(fils,fnams));
        if k > length(avprms)
            avprms(k).corrupt = {};
        end
    catch
        flds = {'subj' 'cond_subj' 'cond' 'tsk' 'alphs' 'bprm1' 'bprm2' 'ls' 'fls' 'iai' 'corrupt'};
        avprms = cell2struct(cell(183*4,length(flds)),flds,2);%183x4 subjects x subconditions
        counter = 0;
    end
else
    outname = [];
    flds = {'subj' 'cond_subj' 'cond' 'tsk' 'tsk' 'alphs' 'bprm1' 'bprm2' 'ls' 'fls' 'iai' 'corrupt'};
    avprms = cell2struct(cell(183*4,length(flds)),flds,2);
    counter = 0;
end



%% go
% for rate = rates
while ~finito%#ok
    [avprms,counter] = Do1Sbj(NAMES, cnd, subj,avprms,counter);
    [cnd, subj, finito] = JacoClock(amnt_sbjcts, cnd, subj);    % advaning us in Jaco clock
end
% end

%% assisting functions


function [avprms,counter] = Do1Sbj(NAMES, cnd, subj, avprms,counter) % NOTE: consider making task_flag an input
global out_paths subconds 

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
     teeg = RemoveHB1Task(i,final,Comps2Reject)'; % using {i} to avoid erasing task seperation, to be differentiated if needed
     [smpls,sensr] = size(teeg);
     epoch_N = smpls/385;
     teeg = reshape(teeg,[385 epoch_N sensr]);
     if z_flag
         teeg = zscore(teeg);
     end
     teeg = cat(1,teeg,zeros([1,epoch_N,sensr]));
     teeg = reshape(teeg,[epoch_N*386,sensr]);
     
     counter = counter + 1;
     avprms(counter).cond = cnd;
     avprms(counter).subj = ceil(counter/length(subconds));
     avprms(counter).cond_subj = subj;
     avprms(counter).tsk = i;
     try
         [avprms(counter).alphs,avprms(counter).bprm1,avprms(counter).bprm2,...
             avprms(counter).ls,avprms(counter).fls,avprms(counter).iai] = eeg2avalnch_epoch(teeg,1:10,[],2,385);
         avprms(counter).corrupt = 0;
     catch
         avprms(counter).corrupt = 1;
     end
     toc % FIXME: missing tic
end



clear final Comps2Reject   % just in case...

function [DATAnhb] = RemoveHB1Task(i,final,Comps2Reject)
global subconds num lim
DATAwhb.data = final.(subconds{i}).data;
DATAwhb.w = final.(subconds{i}).w;
DATAwhb.sph = final.(subconds{i}).sph;
rejcomps = Comps2Reject.(subconds{i});

DATAnhb = remove_HB(DATAwhb, rejcomps, num, lim);

clear DATAwhb rejcomps

% DATAnhb = reshape(DATAnhb,206,385,[]);



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
