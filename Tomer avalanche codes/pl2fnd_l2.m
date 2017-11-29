function [bs,as] = pl2fnd_l2(f,bn,xmin,xmax)

if ~exist('bn','var') || isempty(bn)
    bn = 0;
end

if ~exist('xmax','var') || isempty(xmax)
    if bn
        xmax = find(f,1,'last');
    else
        xmax = round(2/3*find(f,1,'last'));
    end
end

if ~exist('xmin','var') || isempty(xmin)
    xmin = 1;
end
if ~bn
    ys = squeeze(log(f(xmin:xmax)));
    xs = log(xmin:xmax)';
else
    [ys,xs] = hist2smoth(squeeze(f(1:xmax)));
    msk = xs>=xmin;
    ys = log(ys(msk));
    xs = log(xs(msk)');
end
inf_msk = ~isinf(ys);
cc = polyfit(xs(inf_msk),ys(inf_msk),1);
bs = - cc(1);
as = cc(2);