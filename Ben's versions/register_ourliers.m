% Register Outliers
%% replacing nan values with MEANS of params
badsets = []; % preallocating
for i = param_rows
    row = (sts(i,:));
    m = mean(row,'omitnan');
    n = find(isnan(row));
    sts(i,n) = m;
    badsets = [badsets n];
end
badsets = unique(badsets)'; % saving to erase from general analysis (I don't want to add artificial data)

%% looking at extreme sts vals

zz = zscore(sts,1,2); % along rows
ext = zz>= ext_STD | zz<=-ext_STD;

%% Registering
ext_inds = find(ext);
siz = size(ext);
ext_vals = zz(ext);
subsc = inds2subs (siz, ext_inds);

%% in case we are only running analysis on a certain cond, and not all,
% we must align the indeces of the subsc, yo get the right 'badsubj':
addrows = 0;
if cond_flag ~= 0 || cond_flag ~= 1 %(it does not matter if only the first cond is analyzed, or all conds)
    switch cond_flag
        case 2
            addrows = length(NAMES{1})*4; % 4 = # tasks.
        case 3
            addrows = (length(NAMES{1}) + length(NAMES{2}) )*4;
        case 4
            addrows = (length(NAMES{1}) + length(NAMES{2}) + length(NAMES{3}) )*4;
    end
end

subsc(:,2) = subsc(:,2) + addrows;

badsubjs = cell(size(subsc,1),4); % 4 = parameter, name, task, value.
for i = 1:size(subsc,1)
    badsubjs{i,1} = nms{subsc(i,1)};
    [badsubjs{i,2}, badsubjs{i,3}] = watsubj(subsc(i,end));
    badsubjs{i,4} = ext_vals(i);
end
badsubjs = sortrows(badsubjs,2);
unibads = sort(unique({badsubjs{:,2}}'));


% cd(avl_outpath)
% SaveUnique('with outliers')
%% iinput extreme vals!
% [stsrow, matrix] = findinsts(12.69,sts);
% [subj, task] = watsubj(matrix);


%% inner functions

function subs = inds2subs(size,inds)
[s1,s2] = deal(nan(length(inds),1));
for i = 1:length(inds)
    [s1(i),s2(i)] = ind2sub(size,inds(i));
end
subs = cat(2,s1,s2);
end
