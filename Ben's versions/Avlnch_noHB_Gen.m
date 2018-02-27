%LZC_noHB_Gen

clear
clc
start_ben
global out_paths conds subconds num lim plothb E_T_flag  %#ok<NUSED>
DOC_basic
Avlnch_noHB_param

[f,bprm,bprm2,ls,fls,iai] = Avlnch_noHB(data_ratio,task_flag);

%% plotting
