clear all 
%% import 

cd('/Users/admin/Dropbox/Ben Serota/eeg ANALYSES/results/avalanches/important WS/final multis')
% alphas
load('Multi_STS_CTRL__2018_5_28_9_51_18','Alphas','Sigmas');
al_ctr = Alphas;
sig_ctr = Sigmas;

load('Multi_STS_EMC__2018_5_28_12_55_49','Alphas','Sigmas');
al_emcs = Alphas;
sig_emcs = Sigmas;

load('Multi_STS_MC__2018_5_28_12_8_54','Alphas','Sigmas');
al_mcs = Alphas;
sig_mcs = Sigmas;

load('Multi_STS_VS__2018_5_28_10_10_39','Alphas','Sigmas');
al_vs = Alphas;
sig_vs = Sigmas;

clear Alphas Sigmas

%% analyze

ALPHAS =  [mean(al_ctr,1),mean(al_emcs,1),mean(al_mcs,1),mean(al_vs,1)];
SIGMAS =  [mean(sig_ctr,1),mean(sig_emcs,1),mean(sig_mcs,1),mean(sig_vs,1)];
grp = [ones(1,10),2*ones(1,10),3*ones(1,10),4*ones(1,10)];
[al_P,al_t,stats,~] = anovan(ALPHAS,grp','sstype',1,'model','full','varnames',{'condition'},'display','off');
[sig_P,sig_t,stats,~] = anovan(SIGMAS,grp','sstype',1,'model','full','varnames',{'condition'},'display','off');

%prep for alpha cell
ALPHAS_cell =  {mean(al_ctr,1),mean(al_emcs,1),mean(al_mcs,1),mean(al_vs,1)};
SIGMAS_cell =  {mean(sig_ctr,1),mean(sig_emcs,1),mean(sig_mcs,1),mean(sig_vs,1)};

[al_H, al_Pt, al_p_inds] = BensTtest2(ALPHAS_cell,0.05);
[sig_H, sig_Pt, sig_p_inds] = BensTtest2(SIGMAS_cell,0.05);

% ACCESSORY 

function [H,P,inds] = BensTtest2(data, alpha)
% calcs 6 pairs of t-tests and returns p value

[H,P] = deal(nan(length(data)));


for i = 2:4
    [H(1,i), P(1,i)] = ttest2(data{1},data{i},'Vartype','unequal','Alpha',alpha);
end

for i = 3:4
    [H(2,i),P(2,i)] = ttest2(data{2},data{i},'Vartype','unequal','Alpha',alpha);
end

[H(3,4),P(3,4)] = ttest2(data{3},data{4},'Vartype','unequal','Alpha',alpha);

inds = find(~isnan(P));
P = P(inds);
H = H(~isnan(H));
end