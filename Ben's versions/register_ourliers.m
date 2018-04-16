% Register Outliers
ext_inds = find(ext);
siz = size(ext);
ext_vals = zz(ext);
subsc = inds2subs (siz, ext_inds);
% subsc = sortrows(subsc,1);

badsubjs = cell(length(subsc),4); % 4 = parameter, name, task, value.
for i = 1:length(subsc)
    badsubjs{i,1} = nms{subsc(i,1)};
    [badsubjs{i,2}, badsubjs{i,3}] = watsubj(subsc(i,end));
    badsubjs{i,4} = ext_vals(i);
end
badsubjs = sortrows(badsubjs,2);
unibads = unique({badsubjs{:,2}}');
cd(avl_outpath)
SaveUnique('with outliers')
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
