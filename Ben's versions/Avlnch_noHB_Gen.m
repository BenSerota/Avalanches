%LZC_noHB_Gen

% clear
close all
clc
start_ben
global out_paths conds subconds num lim z_flag avl_outpath %#ok<NUSED>
DOC_basic
Avlnch_noHB_param

%% run analysis or load data
if ~exist('avprms','var')  % check if doc_avprms2sts has already run
    cd(avlnch_rslts)
%     try
%         load(outname)
%     catch
        [avprms] = Avlnch_noHB(data_frac, outname);
%     end
end
%% testing resulting avlnchs stats
amnt_crpt = 0;
for i = 1:length(avprms)
    amnt_crpt = amnt_crpt + avprms(i).corrupt;
end
if amnt_crpt>0
    error('some rows are curropt');
end
%% Generate statistics
if ~exist('sts','var')  % check if doc_avprms2sts has already run
    [sts,nms,tsk,grp,subj,hs,phsps] = doc_avprms2sts(avprms,'all',1,5);
end

%% plotting
scatter_prms_by_cond(sts,grp);


%% inner functions
function [] = scatter_prms_by_cond(prms,grp)
global conds
params = {'sigmas', 'alphas', 'taus', 'gammas', 'deltas', 'kappas', 'genkappas'};
param_rows = [2,4,6,8,10,31,32];

for i = 1:length(params)
    figure('name',params{i})
    row = param_rows(i);  % parameter row
    for ii = 1:length(conds) %conds
        toscatter = prms(row , grp == ii);
        
        % throwing away outliers (lower and upper 5prcnt)
%         low = prctile(toscatter,5);
%         up = prctile(toscatter,95);
%         toscatter = toscatter(toscatter>=low & toscatter<=up);
%         
        scatter(ii*ones(1,length(toscatter)),toscatter)
        means(ii) = mean(toscatter);
        hold on
        scatter(ii,means(ii),'bd')
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

% for i = 1:length(conds)
%     sigmas{i} = prms(2,grp==i);
%     alphas{i} = prms(4,grp==i);
%     taus{i} = prms(6,grp==i);
%     gammas{i} = prms(8,grp==i);
%     deltas{i} = prms(10,grp==i);
%     kappas{i} = prms(31,grp==i);
%     genkappas{i} = prms(32,grp==i);
% end

