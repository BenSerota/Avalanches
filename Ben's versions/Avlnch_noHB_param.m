% LZC_noHB_param
%(data_ratio, rates, rate_flag, event_flag ,task_flag);

%% direct inputs
% how much of data to take into consideration (~(0,1))?
data_frac = 1; %0.042; % (rate min = 0.042)

% name of saved file
outname = ['AllSubj_AvlnchPrm' '_thresh3to4'];
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
% 11 thresholds X 10 Time bins.
thresh = 3.1:.1:5; % 2:.1:3;
% thresh = 0.005;
tb_size =  1:10; %1:5; 
pos = 2;


%% avlnch statistics 
% what is a too sparse data set? below __ samples:
sprs = 30000;
% what constitutes an outlier (in STD)?
ext_STD = 4;

%parameters:
params = {'sigmas', 'alphas', 'taus', 'gammas', 'deltas', 'kappas', 'genkappas', 'CutOffs'};

%data in sts matrix in rows:
param_rows = [2,4,6,8,10,31,32,33];
% alpha = sts(4,:);
% tau = sts(6,:);
% sigma = sts(2,:);
% gamma = sts(8,:);
% delta = sts(10,:);
% kappa = sts(31,:);
% kappagen = sts(32,:);

%% testing for significance, at alpha = 
alph = 0.05;

% running multi analysis?
multi_flag = 1;

% multi analysis? over condition # ? (1=VS/2=MC/3=EMC/4=CTRL , 0 = all)
cond_flag = 4; 
