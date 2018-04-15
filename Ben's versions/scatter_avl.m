
%% plotting
bads = prblmdata(30000);
bads(bads>length(sts)) = []; % getting rid of any corrupt data sets. 
sts(:,bads) = [];
grp(:,bads) = [];

%% looking at extreme sts vals
zz = zscore(sts,1,2); % along rows
ext = zz>= ext_STD | zz<=-ext_STD;

%% plot
scatter_prms_by_cond(sts,grp,ext);
tilefigs

%% back to finding extreme subjects names
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
save('New avl from 2 to 3')
%% iinput extreme vals!
% [stsrow, matrix] = findinsts(12.69,sts);
% [subj, task] = watsubj(matrix);


%% inner functions
function [] = scatter_prms_by_cond(prms,grp,ext)
global conds
params = {'sigmas', 'alphas', 'taus', 'gammas', 'deltas', 'kappas', 'genkappas'};
param_rows = [2,4,6,8,10,31,32];

for i = 1:length(params) % over parameters
    figure('name',params{i})
    row = param_rows(i);  % parameter row
    for ii = 1:length(conds) %conds
        toscatter = prms(row , grp == ii);
        scatter(ii*ones(1,length(toscatter)),toscatter)
        means(ii) = mean(toscatter);
        hold on
        scatter(ii,means(ii),'bd')
        extremes = toscatter(ext(row,grp==ii));
        scatter(ii*ones(1,length(extremes)),extremes,'k*')
        
        % setting figure borders
        uplim(ii) = max(toscatter);
        lowlim(ii) = min(toscatter);
    end
    plot(1:4,means)
    title([params{i} ' per lvl of Consciousness'])
    xticks(1:4)
    set(gca,'xticklabel',conds)
    xlabel('Level of Consciousness')
    xlim([0.5 4.5])
    uplim = max(uplim);
    lowlim = min(lowlim);
    ylim([lowlim-.2 uplim+.2]) % 
end
tilefigs
end

function subs = inds2subs(size,inds)
[s1,s2] = deal(nan(length(inds),1));
for i = 1:length(inds)
    [s1(i),s2(i)] = ind2sub(size,inds(i));
end
subs = cat(2,s1,s2);
end
