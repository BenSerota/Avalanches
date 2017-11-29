function [] = loggingplots (SUBJECTS)

% calc alphas
PLOTS.alphas = mean(cell2mat(cellfun(@(x) x.alphas, SUBJECTS, ...
    'UniformOutput', 0)), 1);

Probs = {'Ps',  'Pl',  'Psl', 'Pi'};
for i = 1:length(Probs)
    P = Probs{i};
    PLOTS.(P) = mean(cell2mat(cellfun(@(x) GetProbs(x, P), SUBJECTS, ...
        'UniformOutput', 0)), 1); % mean is according to 1st dimension which is subjects. 
    % eliminating dummy (1st) dimension, after having averaged across
    % subjects:
    PLOTS.(P) = reshape(PLOTS.(P), size(PLOTS.(P), 2), size(PLOTS.(P), 3));
    PLOTS.(P) = (PLOTS.(P))'; % we trandpose because later loglog function
    % will receive a matrix and will plot every row as a different line. yay
end

% creating labels
labels = cellstr(num2str((1:size(SUBJECTS{1,1}.size,2))', 'tb=%-d'));

% Finally plotting
figure('Name','Group Probability of Avalanche by Size, log-log');
loglog(PLOTS.Ps);
title('Group Probability of Avalanche by Size, log-log')
xlabel('Avalanch Size')
ylabel('Probability')
legend(labels)

figure('Name','Group Probability of Avalanche by Length, log-log');
loglog(PLOTS.Pl);
title('Group Probability of Avalanche by Length, log-log')
xlabel('Avalanch Length')
ylabel('Probability')
legend(labels)

figure('Name','Group Probability of Avalanche sizes by Length, log-log');
loglog(PLOTS.Psl)
title('Group Distribution of Avalanche Sizes by Length, log-log')
xlabel('Avalanch Size')
ylabel('Avalanch Length')
legend(labels)

figure('Name','Group Probability of Inter-Avalanch-Interval by Length , log-log');
title('Group Probability of Inter-Avalanch-Interval by Length , log-log')
loglog(PLOTS.Pi)
xlabel('Interval Length')
ylabel('Probability')
legend(labels)

%plotting alphas
figure('Name','Group Alpha values per Time-Bin size');
plot(PLOTS.alphas)
hold on
refline(0,-1.5,'-g') % refeence for wanted alpha value
title('Group Alpha values per Time-Bin size')
xlabel('Time-Bin Size')
ylabel('Alpha value')

tilefigs

function P = GetProbs(BIG5, field)

switch field
    case 'Ps'
        f = 'size';         % for avalanch SIZES
    case 'Pl'
        f = 'length';       % for avalanch LENGTH
    case 'Psl'
        f = 'size_length';	% for avalanch SIZES per LENGTH
    case 'Pi'
        f = 'iai_length';	% for IAI
    otherwise
        error('Unknown field, bitch');
end

N = size(BIG5.(f),2);       % number of time-bin sizes
K = length(BIG5.(f){1,1});	% TODO: find out what this sparse format is / is good for
P = zeros(1, N, K);

for i = 1:N % looping over tb sizes
    A = cat(2,BIG5.(f){:,i});
    temp_i = sum(A,2);
    P(1,i,:) = full(temp_i./sum(temp_i));
end
