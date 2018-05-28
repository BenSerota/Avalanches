% PlotLogs

load('NotAve_CTRL_thresh_3_4p5.mat')
% load('Multi_STS_VS__2018_5_28_10_10_39.mat')

D4Plots = cellfun(@(x) mean(x,1),avprms,'uniformoutput',0);
A = reshape([avprms.alphs],[],10);
A = cellfun(@(x) mean(x,1),A,'uniformoutput',0);

t = avprms(exmp).ls{exmp, exmp}';
T = sum(t);
loglog(1:length(t),t/T)