function [BIG5] = eeg2avalnch_ben(raw_data,tb_size,thresh,pos,cont)
%function [av_size,av_sigma,bprm2,av_length,av_size_length,iai_length,wv] = eeg2avalnch_ben(raw_data,tb_size,thresh,pos,cont)

%   derive avalanches from multi electrode EEG time series (column)
%   eeg2avalnch3(EEG,tbin,thresh,{0|1|2} for neg/pos/both,control_flag - 1 for control distrib);
%   INPUTS:
%    ts = data time x sensors (elec = COLUMN)
%    dt = time bin size (i.e. size of window that checks if event has occured,
%           always returning 1|0) e.g. [1:5].
%   thresh = 2 options either STD or in [] which takes the 99th prcntl +-
%   0.2. STD, in incraments of 0.1., and averages across. because 99th
%   prcntl=> about 2.8 STD, the threshholds will be (in std): 
%   -3, -2.9, -2.8, -2.7, -2.6 (average) and 2.6, 2.7, 2.8, 2.9, 3 (average). 
%   
%   pos = 0|1|2 . 0 = neg, 1 = positive, 2 = both.
%   cont = flag for control. 0|1. 0 = no, 1 = yes. control shows exp dist
%   is not trivial.

% OUTPUTS:
% f = avlnch size (hist of sum between two 0s) (histogramable)
% bprm = sigma
% bprm2 = a new sigma calc way, for undersampling. not used right now.
% ls = avlnch duration dist (histogramable)
% fls = dist of avlnch size per avlnch length (duration)
% iai = inter avalanch interval length (dist of quiet times, like truce times)
% wv = shape of avalanche,

% TWIST OF FAITH:
% there is a possible offset to every time-bin(tb) size: none for tb=1,
% 1 for tb=2 etc... i.e. for tb=1 there is 1 option, for tb=2 there are 2 etc. 
% therefore, for dt=3, we have 1+2+3 manifestations, so we will implicitly 
% have 6 f values (and histograms in general), over which we average!.
% in the same sense, for dt = 1:5 we have 15 implicit manifestations.

% 


bprm2 = nan(length(thresh),length(tb_size));
[av_size,frst,scnd,av_length,av_size_length,iai_length] = deal(cell(length(thresh),length(tb_size)));
if nargout == 7
    wv = cell(length(thresh),length(tb_size));
end
% a flag for computing a control distribution through permutation
if ~exist('cont','var')
    cont = 0;
end
if isempty(thresh)
    thresh = [zeros(1,5);-0.2:0.1:0.2];
end

data = zscore(raw_data);

for t = 1 : size(thresh,2)
    msk = ts2ev(data,thresh(:,t),pos);
    for d = 1 : length(tb_size)
        if nargout == 8
            [av_size{t,d},bprm2(t,d),frst{t,d},scnd{t,d},av_length{t,d},av_size_length{t,d},iai_length{t,d},wv{t,d}] = msk2f(msk,tb_size(d),cont);
        else
            [av_size{t,d},bprm2(t,d),frst{t,d},scnd{t,d},av_length{t,d},av_size_length{t,d},iai_length{t,d}] = msk2f(msk,tb_size(d),cont);
        end
    end
end

av_sigma = cellfun(@(x,y) mean(y./x),frst,scnd);

% BEN ADDED
alphas = cellfun(@(x) pl2find_ML(x),av_size); %TODO: is this right?? 
% A. is this missing a minus? values tend to be between .5-2 .
% B.should I not calculate this regarding each form of avalanch? 
%i.e. for each size/length/iai/sizeXlength ?

%into a structure
BIG5.alphas = mean(alphas,1);
BIG5.size = av_size;
BIG5.length = av_length;
BIG5.size_length = av_size_length;
BIG5.iai_length = iai_length;


function msk = ts2ev(ts,thresh,pos)
%Assumes columns
if ~thresh(1)
    if length(thresh) == 1
       thresh(2) = 0; 
    end
    dstrb = sort(ts(:));
    l = length(dstrb);
    if ~pos
        k = round(l*0.01);
        thresh(1) = dstrb(k);
        msk = ts < (thresh(1)-thresh(2));
    elseif pos == 1
        k = round(l*0.99);
        thresh(1) = dstrb(k);
        msk = ts > sum(thresh);
    elseif pos == 2
        k = round(l*0.995);
        thresh(1) = dstrb(k); % print this 
        msk = ts > sum(thresh);
        k2 = round(l*0.005);
        thresh2 = dstrb(k2); % print this 
        msk = msk | (ts < (thresh2-thresh(2)));
    end
else
    if ~pos
        msk = ts < -thresh;
    elseif pos == 1
        msk = ts > thresh;
    elseif pos == 2
        msk = abs(ts)>thresh;
    end
end
smpls = size(msk,1);
msk = col([msk;zeros(1,size(ts,2))]);
ts = col([ts;zeros(1,size(ts,2))]);
[s,e] = enpoints2find(msk);
if isempty(s)
    return
end
mx = max(e-s+1);
inds = round(linspacen(s,e,mx));
vals = ts(inds);
if ~pos
    [~,inds2] = min(vals);
elseif pos == 1
    [~,inds2] = max(vals);
elseif pos == 2
    [~,inds2] = max(abs(vals));
end
inds2 = sub2ind(size(vals),inds2,1:size(vals,2));
inds = inds(inds2);
msk = msk*0;
msk(inds) = 1;
msk = reshape(msk,[smpls+1 length(msk)/(smpls+1)]);
msk(end,:) = '';

function [f2,bprm,fst,scnd,ls,fls,iai,wv] = msk2f(msk,dt,cont,rf)
%avlnch size counts,auto-regressive branching parameter estimate, 1st event sizes,2nd event size,duration counts,sum of
%events sizes per duration,IAI, wave form

%rf is the recursion flag
if ~exist('rf','var')
    rf = 1;
end
if dt > 1 && rf
    [f2,bprm,ls,fls,iai,wv] = deal(sparse(1,1));
    [fst,scnd] = deal([]);
    for k = 1 : dt
        if nargout == 8
            [tf2,tbprm,tfst,tscnd,tls,tfls,tiai,twv] = msk2f(msk(k:end,:),dt,cont,0);
        else
            [tf2,tbprm,tfst,tscnd,tls,tfls,tiai] = msk2f(msk(k:end,:),dt,cont,0);
        end
        f2 = f2 + (1/dt)*tf2;
        bprm = bprm + (1/dt)*tbprm;
        fst = [fst; tfst];
        scnd = [scnd; tscnd];
        iai = iai + (1/dt)*tiai;
        ls = ls + (1/dt)*tls;
        fls = fls + (1/dt)*tfls;
        if nargout == 8
            wv = wv + (1/dt)*twv;
        end
    end
    return
end
%bin by steps of dt and leave up to one event only per sensor in bin
[smpls,elctrd] = size(msk);
N = floor(smpls/dt);
if dt > 1
    msk = reshape(msk(1:N*dt,:),[dt N elctrd]);
    msk = squeeze(sum(msk))>0;
end
if cont
    [~,perms] = sort(rand(size(msk)));
    perms = bsxfun(@plus,perms,(0:elctrd-1)*N);
    msk = msk(perms);
end
f = sum(msk,2);
bprm = sigma_MR_estimator(row(f));
f2 = sparse(35000,1);
ls = sparse(5000,1);
fls = sparse(5000,1);
iai = sparse(5000,1);
if nargout == 8
    wv = zeros(100);
end

ee = find(~f);
if isempty(ee)
    [fst,scnd] = deal(nan);
    return
end
%truncate avalanches that weren't captured in full
%at the beginning
if ee(1)~=1
    f(1:ee(1)-1) = 0;
end
%and at the end
if ee(end)~= N
    f(ee(end) : N) = 0;
end

[s,e] = enpoints2find(f>0);

%find the 1st and 2nd event size for each avalanche
fst = f(s);
scnd = f(s+1);

%compute event size, duration, size per duration
msk = cumsum(f);

%avln durations
tls = e-s+1;

%avlnch magnitudes
es = diff([0;msk(e)]);
counts = diff([0;find(diff(sort(es))~=0);length(es)]);
f2(unique(es)) = counts;
[tls,ord] = sort(tls);
[tls,durs] = unique(tls,'last');
ls(tls) = diff([0;durs]);

%avlnch mags per durs
es = cumsum(es(ord));
fls(tls) = diff([0;es(durs)]);

%IAI
tiai = s(2:end)-e(1:end-1);
[tiai,ord] = sort(tiai);
[tiai,ais] = unique(tiai,'last');
iai(tiai) = diff([0;ais]);

%find avlnch waveforms
if nargout == 8 && ~isempty(s)
    twv = interp1((1:length(f))',f,linspacen(s,e,100));
    twv = cumsum(twv(:,ord),2);
    wv(:,tls) = diff([zeros(100,1) twv(:,durs)],[],2);
end

function [s,e] = enpoints2find(msk)
%msk is assumed to be a column binary vector
dm = diff(msk);
e = find(dm==-1);
s = find(dm==1)+1;
if msk(1)
    s = [1;s];
end
if msk(end)
    e(end+1) = length(msk);
end

function x = col(x)

x = x(:);

function x = row(x)

x = x(:)';

