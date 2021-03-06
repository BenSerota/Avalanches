%% replacing nan values with MEANS of params

% assembling defects from 3 origins :
% 1. failed avalanche analysis
% 2. too sparse
% 3. crazy recording fuck up (discovered via outlier)

% 1. failed avalanche analysis
%  = badsets

        % badsets = []; % preallocating
        % for i = param_rows
        %     row = (sts(i,:));
        %     m = mean(row,'omitnan');
        %     n = find(isnan(row));
        %     sts(i,n) = m;
        %     badsets = [badsets n];
        % end
        % badsets = unique(badsets)'; % saving to erase from general analysis (I don't want to add artificial data)

% 2. too sparse
spdata = sparsedata(sprs); % bads = too sparse data sets

% 3. crazy recording fuck up (discovered via outlier)
% restoring subsc to original size of specific matrix
spec_subsc = subsc;
spec_subsc = spec_subsc(:,2) - addrows;
f_ups = unique(spec_subsc); % = f up (recording?) data sets

% assembling
bads = [badsets ; spdata ; f_ups]; % adding sets which failed avalanche analysis (for some reason)
bads = sort(bads);
bads = unique(bads);


% erasing all
avprms(bads) = [];
sts(:,bads) = [];
grp(:,bads) = [];
alphs(:,bads) = [];
sigms(:,bads) = [];
taus(:,bads) = [];
gams(:,bads) = [];
delts(:,bads) = [];
kap(:,bads) = [];
gen_kp(:,bads) = [];
cos(:,bads) = [];