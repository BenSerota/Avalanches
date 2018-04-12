function [name, task] = watsubj(num)
start_ben; 
DOC_basic;
DOC_basic2;
all = 732;
vs = 77;
mcs = 70;
emcs = 24;
ctrl = 12;

subj = floor(num/4);
task_n = mod(num,4);
if task_n == 0
    task_n = 4;
end
task = subconds{task_n};

if subj <= vs % subj is in VS
    name = NAMES{1}{subj};
elseif subj > vs && subj <= vs+mcs % subj is in MCS
    subj = subj - vs;
    name = NAMES{2}{subj};
elseif subj > vs+mcs && subj <= vs+mcs+emcs % subj is in EMCS
    subj = subj - (vs+mcs);
    name = NAMES{3}{subj};
elseif subj > vs+mcs+emcs && subj <= all % subj is in CTRL
    subj = subj - (vs+mcs+emcs);
    name = NAMES{4}{subj};
else
    error('number is out of range')
end

