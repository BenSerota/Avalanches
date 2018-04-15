% LZC_noHB_param
%(data_ratio, rates, rate_flag, event_flag ,task_flag);

%% direct inputs
% how much of data to take into consideration (~(0,1))?
data_frac = 1; %0.042; % (rate min = 0.042)

% name of saved file
outname = ['AllSubj_AvlnchPrm' '_thresh2to3'];
avl_outpath = 'C:\Users\BenSerota\Dropbox\Ben Serota\eeg ANALYSES\results\avalanches';

%% globals
%Z scoring. 1 = per epoch, 0 = per row, concatenated all epochs
z_flag = 1;

% r chooses how many of the first ICs to reject in each dataset:
num = 30;

% what's the component beyond which we don't care if there is HB (as
% components are organized hierarchically).
lim = 30;

%% avalanch [arameter values
tb_size = 1:10; 
thresh = 2:.1:3;
% thresh = 0.005;
pos = 2;
% cont = 1;

%% avlnch statistics 
% what constitutes an outlier (in STD)?
ext_STD = 4;
