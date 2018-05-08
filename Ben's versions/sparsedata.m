function [vec] = sparsedata(cutoff)
global cond_flag
DOC_basic
cd(LZ_rslts)
load('numelements')
if cond_flag == 1 || cond_flag == 2 || cond_flag == 3 || cond_flag == 4
    elements = {elements{cond_flag}};
end

elements = cellfun(@(x) x', elements,'uniformoutput', false);
elements = cellfun(@(x) x(:), elements,'uniformoutput', false);
% now data is 1234 1234 1234 etc... in terms of task.
allmat = cat(1, elements{:});
badinds = bsxfun(@lt,allmat,cutoff);
vec = find(badinds);