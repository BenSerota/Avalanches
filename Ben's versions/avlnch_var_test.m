%avlnch var test
prms = {'alphas', 'taus', 'sigma', 'gamma', 'delta', 'kappa', 'kappagen'};

% alpha = sts(4,:);
% tau = sts(6,:);
% sigma = sts(2,:);
% gamma = sts(8,:);
% delta = sts(10,:);
% kappa = sts(31,:);
% kappagen = sts(32,:);

codes = [4,6,2,8,10,31,32];

for i = 1:length(prms)
    for ii = 1:length(conds)
        prms{i}.(conds{ii}) = sts(codes(i,[avprms.cond]==ii));
    end
end