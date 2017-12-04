function kappa = kappa_est(dist,bin)
%KAPPA_EST Estimate Kappa value for a cascade size distribution
%use the distribution binned to 10 logarithmic bins
% compute empirical CDF

if ~exist('bin','var')
    bin = 0;
end

if size(dist,2) > 1
    n = size(dist,2);
    kappa = nan(n,1);
    for k = 1 : n
        kappa(k) = kappa_est(dist(:,k),bin);
    end
    return
end

if bin
    [dist,~,bin_edgs] = hist2smoth(dist,find(dist,1,'last'),20);
else
    l = find(dist,1,'last');
    bin_edgs = (.5:l+.5);
    dist = dist(1:l);
end
bin_width = diff(bin_edgs');
F_emp = cumsum(dist.*bin_width); % CDF
F_emp = F_emp./F_emp(end);

% compute theoretical CDF
l = bin_edgs(1);
L = bin_edgs(end);
F_theo1 = (1 - sqrt(l./bin_edgs')) / (1 - sqrt(l/L));
F_theo = F_theo1(2:end);

kappa = 1 + mean(F_theo - F_emp);


end

