function [vec] = prblmdata(cutoff)
DOC_basic
cd(LZ_rslts)
load('numelements')
elements = cellfun(@(x) x', elements,'uniformoutput', false);
elements = cellfun(@(x) x(:), elements,'uniformoutput', false);
% now data is 1234 1234 1234 etc... in terms of task.
allmat = cat(1, elements{:});
badinds = bsxfun(@lt,allmat,cutoff);
vec = find(badinds);