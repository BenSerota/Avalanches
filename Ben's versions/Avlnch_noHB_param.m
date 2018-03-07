% LZC_noHB_param
%(data_ratio, rates, rate_flag, event_flag ,task_flag);

%% direct inputs
% how much of data to take into consideration (~(0,1))?
data_frac = 1; %0.042; % (rate min = 0.042)

%use zscore only(=0) or hilbert transform on zscore(=1)?
event_flag = 1; 

% should we also seperate analysis by task (0=no)?
task_flag = 0;


%% globals
% r chooses how many of the first ICs to reject in each dataset:
num = 1;

% what's the component beyond which we don't care if there is HB (as
% components are organized hierarchically).
lim = inf;

% calc LZC by row (temporal; 0) or column (sptial; 1)
E_T_flag = 0;

% % plot HB EEG ?
% plothb = 0;
