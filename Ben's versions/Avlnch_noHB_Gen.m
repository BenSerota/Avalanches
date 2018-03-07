%LZC_noHB_Gen

clear
clc
start_ben
global out_paths conds subconds num lim plothb E_T_flag  %#ok<NUSED>
DOC_basic
Avlnch_noHB_param

[ALL_BIG5s] = Avlnch_noHB(data_frac,task_flag);

%% plotting
