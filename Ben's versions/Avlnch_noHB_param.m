% LZC_noHB_param
%(data_ratio, rates, rate_flag, event_flag ,task_flag);

%% direct inputs
% how much of data to take into consideration (~(0,1))?
data_frac = 1; %0.042; % (rate min = 0.042)

% path of saved file is 'avlnch_rslts' in DOC_basic;

%% globals
%Z scoring. 1 = per epoch, 0 = per row, concatenated all epochs
z_flag = true;

% r chooses how many of the first ICs to reject in each dataset:
num = 30;

% what's the component beyond which we don't care if there is HB (as
% components are organized hierarchically).
lim = 30;

%% avalanch [arameter values
%  thresholds X Time bins.
% thresh =  2:.1:3; %3:.1:4; %.1:.1:5; 0.005;
thresh = 3:.1:4.5;
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
multi_flag = 0;

if ~multi_flag % for single thresh and tb statistics:
    chosen_th = 5; % chosen threshold index (out of thresh)
    chosen_tb = 1; % chosen tb (out of tb_size)
end
% multi analysis? over condition # ? (1=VS/2=MC/3=EMC/4=CTRL , 0 = all)
cond_flag = 0;

%name of saved file
if multi_flag
    out_a = 'NotAve_';
else
    out_a = 'Single_';
end

switch cond_flag
    case 1
        out_b = 'VS_';
    case 2
        out_b = 'MC_';
    case 3
        out_b = 'EMC_';
    case 4
        out_b = 'CTRL_';
    case 0
        out_b = 'ALL_CONDS_';
end

out_c = sprintf('thresh_%g_%g',thresh(1),thresh(end));
% if thresholds contain periods: problem for file name, replace with 'p':
out_c = strrep(out_c,'.','p');

outname = [out_a out_b out_c];

LZC_flag = 0; % for statistics, using BensAnovaTest. i.e. doing avalanches not LZC.
