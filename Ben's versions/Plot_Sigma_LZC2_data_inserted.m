
m = [.6 .71; .63 .71; .7 .71; .78 1.05];

plot(m(:,1),m(:,2),'-o','linewidth',2)
title(sprintf('Sigma X LZ Complexity \n Estimation of Integration X Differentiation'))
xlabel('LZC score (differentiation)')
ylabel('Sigma (integration)')
xlim([0.59 .8])
ylim([.7 1.1])
text(m(1,1),m(1,2),'VS','fontsize',14,'color','green', 'fontweight', 'bold') 
text(m(2,1),m(2,2),'MCS','fontsize',14,'color','red', 'fontweight', 'bold') 
text(m(3,1),m(3,2),'EMCS','fontsize',14,'color','magenta', 'fontweight', 'bold') 
text(m(4,1),m(4,2),'CTR','fontsize',14,'color','blue', 'fontweight', 'bold') 
