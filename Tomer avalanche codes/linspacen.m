function y = linspacen(d1, d2, n)
%LINSPACEN Linearly spaced vector.
%   LINSPACEN(X1, X2) generates a row vector of 100 linearly
%   equally spaced points between X1 and X2.
%
%   LINSPACEN(X1, X2, N) generates N points between X1 and X2.
%   For N < 2, LINSPACE returns X2.
%
%   X1, X2 can also be column vectors in which case the output will 
%   be a matrix of spaced columns

if nargin == 2
    n = 100;
end

l1 = length(d1);
l2 = length(d2);
if l1 < l2
    d2(l1+1:end) = '';
elseif l2 < l1
    d1(l2+1:end) = '';
end

spacer = (0:n-2);
step = (d2-d1)/(n-1);
y = [d1(:,ones(1,n-1))+step*spacer d2]';