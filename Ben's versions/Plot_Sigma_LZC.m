load('allLZCs.mat')
load('th=3p4_tb=1.mat')

LZC_vec = nan(732,1); j = 1;
for i = 1:length(LZCs_per_cond)
    h = numel(LZCs_per_cond{i});
    LZC_vec(j:j+h-1) = LZCs_per_cond{i}(:);
    j = j+h;
end

LZC_vec(bads) = [];

bads_n_cond = bads;

%borders between LCs
b1 = 308;
b2 = 588;
b3 = 684;

for i = 1:numel(bads_n_cond)
    if bads_n_cond(i,1) <= b1
        bads_n_cond(i,2) = 1;
    elseif bads_n_cond(i,1) <= b2 
        bads_n_cond(i,2) = 2;
    elseif bads_n_cond(i,1) <= b3
        bads_n_cond(i,2) = 3;
    else 
        bads_n_cond(i,2) = 4;
    end
end



borders = find(diff(bads_n_cond(:,2)));

% new borders, on truncated data:
bb1 = b1- borders(1);
bb2 = b2 - borders(2);
bb3 = b3 - 104;
e = numel(LZC_vec);

% sequences = 
sq1 = 1:bb1;
sq2 = bb1+1:bb2;
sq3 = bb2+1:bb3;
sq4 = bb3+1:e;

scatter(LZC_vec(sq1),sigms(sq1),'g')
hold on
scatter(LZC_vec(sq2),sigms(sq2),'r')
scatter(LZC_vec(sq3),sigms(sq3),'m')
scatter(LZC_vec(sq4),sigms(sq4),'b')

m(1,1) = mean(LZC_vec(sq1));
m(1,2) = mean(sigms(sq1));
m(2,1) = mean(LZC_vec(sq2));
m(2,2) = mean(sigms(sq2));
m(3,1) = mean(LZC_vec(sq3));
m(3,2) = mean(sigms(sq3));
m(4,1) = mean(LZC_vec(sq4));
m(4,2) = mean(sigms(sq4));

plot(m(:,1),m(:,2),'-o','linewidth',2)
title(sprintf('Sigma X LZ complexity \n Estimation of Integration X Differentiation'))
xlabel('LZC score (differentiation)')
ylabel('Sigma (integration)')
text(m(1,1),m(1,2),'VS','fontsize',14,'color','green', 'fontweight', 'bold') 
text(m(2,1),m(2,2),'MCS','fontsize',14,'color','red', 'fontweight', 'bold') 
text(m(3,1),m(3,2),'EMCS','fontsize',14,'color','magenta', 'fontweight', 'bold') 
text(m(4,1),m(4,2),'CTR','fontsize',14,'color','blue', 'fontweight', 'bold') 
