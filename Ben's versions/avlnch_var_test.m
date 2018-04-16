%avlnch var test
prms = {'alphas', 'taus', 'sigma', 'gamma', 'delta', 'kappa', 'kappagen'};

codes = [4,6,2,8,10,31,32];

a = struct();
for i = 1:length(prms)
    for ii = 1:length(conds)
        a = setfield(a,prms{i},conds{ii},...
            sts(codes(i),[avprms.cond]==ii)');
    end
end

% 1. run F test, per parameter
for i = 1:length(prms)
    temp = cell(1,4);
    for ii = 1:4
        temp{ii} = a.(prms{i}).(conds{ii});
    end
    F.(prms{i}) = BensAnovaTest(temp,0.05);
end
% % 2. run paired t-tests
% %     if P <= alpha
% [H, Pt, p_inds] = BensTtest(data,alpha);
% 
% % 3. correct for mult comp
% [Pt_cor, crit_p, h] = fdr_bh(Pt,alpha);
% 
% % 4. prepare P values for Bar graph
% Ps4bar = prepP(Pt_cor,p_inds);
% 
% LZCs_to_bar = cellfun(@(x) mean(x),data,'UniformOutpu',false');
% LZCs_to_bar = cell2mat(LZCs_to_bar);
% save_bar = 0;
% %     E = cellfun(@(x) mean(x,2), STEs);
% E = cell2mat(STEs);
% 
% BensSuperbar(LZCs_to_bar,i,Ps4bar,E,save_bar,LZC_nohb_outpath )
% 
