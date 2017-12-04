function kappa = kappa_gen(dist,bn,meth)
%KAPPA_GEN Estimate Kappa value for a cascade size distribution given a
%fitted exponent

if ~exist('bn','var')
    bn = 0;
end
if ~exist('meth','var')
    meth = 0;
end

if size(dist,2) > 1
    n = size(dist,2);
    kappa = nan(n,1);
    for k = 1 : n
        kappa(k) = kappa_gen(dist(:,k),bn,meth);
    end
    return
end
%find best fit exponent
if meth 
    alph = -pl2fnd_l2(dist,bn);
else
    alph = -pl2find_ML(dist,bn);
end

if bn 
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
xs = (bin_edgs(1:end-1)+bin_edgs(2:end))/2;
F_theo = cumsum(col(xs).^alph);
F_theo = F_theo/F_theo(end);

kappa = 1 + mean(F_theo - F_emp);



