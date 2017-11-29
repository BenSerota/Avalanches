function b = prm2dnoiz(bs,k)
%fit a line to a sequence of crit params (evenly spaced in scale) and
%return the fit for the desired scale
cc = polyfit((1:length(bs))',bs(:),1);
b = polyval(cc,k);