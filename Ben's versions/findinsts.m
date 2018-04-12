function [s1, s2] = findinsts (stsval,sts)
num_dig = digits_debugged(stsval);

linearinds = find(round(sts,num_dig)==stsval);
for i = 1:length(linearinds)
    [s1(i), s2(i)] = ind2sub(size(sts),linearinds(i));
end
% subs = cat(2,s1,s2);

function y = digits_debugged(x)
x = abs(x); %in case of negative numbers
y = 0;
while (floor(x)~=x)
    y = y+1;
    x = x*10;
end