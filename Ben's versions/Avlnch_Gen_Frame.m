% Avlnch_Gen_Frame
% size of av_size/length/size_length histograms is: M x N, where:
%   M = 5 = # of manifestations of thresh, 
%   N = # of different bin sizes

%% for a single subject
% clear 
% close all 
% clc
% 
% av_param_values % sets avalanche parameters
% load setclean %RS1_no_time
% raw_data = t'; %rsEEGec';
% [BIG5] = eeg2avalnch_ben(raw_data,tb_size,thresh,pos,cont);
% loggingplots (BIG5)

%% for multiple subjects ( mult data matrices):

% 
% % set path first
clear 
clc
close all

av_param_values % sets avalanche parameters
fold = what;
subjects = length(fold.mat);
SUBJECTS = cell(subjects,1);

% TODO: maybe size by length is actually only prob of size

for i = 1:subjects
    temp = load(fold.mat{i});
    %strct_names = fieldnames(temp);
    %mtrix = strct_names{1};
    raw_data{i} = temp.rsEEGec'; %temp.t'; % must have rsEEGec structure!
    clear temp
    [BIG4] = eeg2avalnch_ben(raw_data{i},tb_size,thresh,pos,cont); 
    SUBJECTS{i} = BIG4; % feeding into main cell array
    clear BIG4
end

loggingplots(SUBJECTS)

calc_Bentropy(SUBJECTS)