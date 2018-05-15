function avprms = Avlnch_noHB(data_frac,outname)
% this function:
% 1. takes INPUT ratio of jaco data and removes its HB component.
% 2. calculates Avalanches parameters.
% 3. saves output as it goes. able to continue from point at which it stopped.

global out_paths subconds num lim tb_size thresh pos z_flag avlnch_rslts %#ok<NUSED>

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
while ~finito
    [avprms,counter] = Do1Sbj(NAMES, cnd, subj,avprms,counter, z_flag);
    cd(avlnch_rslts)
    save(outname)
    [cnd, subj, finito] = JacoClock(amnt_sbjcts, cnd, subj);                % advaning us in Jaco clock
end



%% assisting functions

function [avprms,counter] = Do1Sbj(NAMES, cnd, subj, avprms,counter, z_flag)
global out_paths subconds tb_size thresh pos

cd(out_paths{cnd})
name = char(NAMES{cnd}(subj));
name_p = [ name '_prep'];
name_I = [ name '_HBICs'];

% load each set of ICs per subj
load(name_p);                                                              % load subject
load(name_I);                                                              % load HB info
tic
for i  = 1:length(subconds)
    teeg = RemoveHB1Task(i,final,Comps2Reject)';
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
            avprms(counter).ls,avprms(counter).fls,avprms(counter).iai] = eeg2avalnch_epoch(teeg,tb_size,thresh,pos,385);
        avprms(counter).corrupt = 0;
    catch
        avprms(counter).corrupt = 1;
    end
end
toc
fprintf('finished subj %s \n', name);

clear final Comps2Reject   % just in case...

function [DATAnhb] = RemoveHB1Task(i,final,Comps2Reject)
global subconds num lim
DATAwhb.data = final.(subconds{i}).data;
DATAwhb.w = final.(subconds{i}).w;
DATAwhb.sph = final.(subconds{i}).sph;
rejcomps = Comps2Reject.(subconds{i});

DATAnhb = remove_HB(DATAwhb, rejcomps, num, lim);

clear DATAwhb rejcomps % just in case...


function [DATAnhb] = remove_HB(DATAwhb, rejcomps, num, lim)
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
