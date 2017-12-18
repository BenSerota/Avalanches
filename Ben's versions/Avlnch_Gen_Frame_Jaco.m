% Avlnch_Gen_Frame
% size of av_size/length/size_length histograms is: M x N, where:
%   M = 5 = # of manifestations of thresh,
%   N = # of different bin sizes

clear
clc
close all
%% Data Analysis

%choose MAC or CPL comp:
%~~~~~~~~~~~
data_paths = {'E:\DOC\Data from Jaco\VS';'E:\DOC\Data from Jaco\MCS';...
    'E:\DOC\Data from Jaco\EMCS';'E:\DOC\Data from Jaco\CTRL'};

mac_data_paths = {'/Users/admin/Desktop/DACOBO_h/VS';'/Users/admin/Desktop/DACOBO_h/EMCS';...
    '/Users/admin/Desktop/DACOBO_h/MCS';'/Users/admin/Desktop/DACOBO_h/CTRL'};
%~~~~~~~~~~~

conds = {'VS';'MCS';'EMCS';'CTRL'};
sbcndt = {'LSGS';'LSGD';'LDGS';'LDGD'};

[tb_size thresh pos cont] = av_param_values; % sets avalanche parameters

for i = 1:length(conds)                                                     % over conditions
    
%     cd(data_paths{i})
    cd(mac_data_paths{i})
    
    fold = what;
    subjects = length(fold.mat);
    SUBJECTS = cell(subjects,1);
    
    % TODO: maybe size by length is actually only prob of size
    
    for ii = 1:subjects                                                     % over subjects
%         temp = load(fold.mat{ii});
        load(fold.mat{ii});
        temp_raw = matcat(temp);
        raw_data{ii} = temp_raw'; % must have correct structure!
        clear temp
        [BIG5] = eeg2avalnch_ben(raw_data{ii},tb_size,thresh,pos,cont);
        SUBJECTS{ii} = BIG5; % feeding into main cell array
        clear BIG4
    end
    
end
%%

loggingplots(SUBJECTS)

calc_Entropy(SUBJECTS)


% % functions
