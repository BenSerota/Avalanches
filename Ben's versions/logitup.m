function [] = logitup(AVPRMS)
% Plots loglog plots of 5 parameters.
params = {'alphs','ls','fls','iai'};
names = {'alphs is = ?', 'ls = LENGTH?', 'fls = SIZE?','Inter-Avalanche-Interval'};
len = length(params);
for i = 1:len
    t = AVPRMS.(params{i}) / sum(AVPRMS.(params{i}));
    len_s = length(AVPRMS.(params{i}));
    series = 1:len_s;
    expseries = series.^-1.5;
    figure()
    loglog(series,t)
    hold on
    plot(linspace(1,len_s,len_s),expseries,'color','red')
    title(sprintf('Probability of Avalanche %s' ,names{i}),'fontsize',14,'FontWeight','bold')
    xlabel(['Avalanche' names{i}], 'fontsize',14,'FontWeight','bold')
    ylabel('Probability','fontsize',14,'FontWeight','bold')
end
tilefigs