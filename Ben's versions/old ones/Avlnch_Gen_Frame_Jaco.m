% Avlnch_Gen_Frame
% size of av_size/length/size_length histograms is: M x N, where:
%   M = 5 = # of manifestations of thresh,
%   N = # of different bin sizes

clear
clc
close all
%% Data Analysis

DOC_basic

[tb_size thresh pos cont] = av_param_values; % sets avalanche parameters
SUBJECTS = cell(1,length(conds)); % just for now. should be subjsXconds.
for i = 1:length(conds)                                                     % over conditions
    
    cd(data_paths{i})
    
    fold = what;
%     subjects = length(fold.mat);
    
    % TODO: maybe size by length is actually only prob of size
    
    for ii = 1:1 %subjects                                                     % over subjects
%         temp = load(fold.mat{ii});
        load(fold.mat{ii});
        temp = data;
        temp_raw = matcat(temp);
        raw_data{ii} = temp_raw'; % must have correct structure!
        clear temp 
        BIG5 = eeg2avalnch_ben(raw_data{ii},tb_size,thresh,pos,cont);
        SUBJECTS{ii,i} = BIG5; % feeding into main cell array
        clear BIG5
    end

end
%%

loggingplots(SUBJECTS)

calc_Entropy(SUBJECTS)
