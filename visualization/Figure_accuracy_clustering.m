%% No Negatives

load('Results_noNegatives.mat');

y1= CC_LFRb - mean(CC_FCsim_abs,2);
y2 = CC_LFRb - mean(CC_FC_perc,2);
y3 = CC_LFRb - mean(CC_FC_objL,2);
y4 = CC_LFRb - mean(CC_FC_stat,2);
y5 = CC_LFRb - mean(CC_FC_tree,2); 

figure; hold on;
plot(y1,'Color',[1 .8 .8],'LineWidth',5) % r for red
plot(y2,'Color',[.6, .6, .6],'LineWidth',4) % r for red
plot(y3,'Color',[1 .4 .4],'LineWidth',3) % r for red
plot(y4,'Color',[.2 .2 .2],'LineWidth',2) % r for red
plot(y5,'Color',[ .4 .4 1 ],'LineWidth',1) % r for red
lgd=legend(NMI_sim_header(2:end));
xlabel('node'); ylabel('LFR CC - group average binary CC');
title('Accuracy of CC Value after Thresholding (No Negatives)');
lgd.FontSize = 16;
ax = gca; 
ax.FontSize = 16; 


%% Absolute Value

load('Results_absValue.mat');

y1= CC_LFRb - mean(CC_FCsim_abs,2);
y2 = CC_LFRb - mean(CC_FC_perc,2);
y3 = CC_LFRb - mean(CC_FC_objL,2);
y4 = CC_LFRb - mean(CC_FC_stat,2);
y5 = CC_LFRb - mean(CC_FC_tree,2); 

figure; hold on;
plot(y1,'Color',[1 .8 .8],'LineWidth',5) % r for red
plot(y2,'Color',[.6, .6, .6],'LineWidth',4) % r for red
plot(y3,'Color',[1 .4 .4],'LineWidth',3) % r for red
plot(y4,'Color',[.2 .2 .2],'LineWidth',2) % r for red
plot(y5,'Color',[ .4 .4 1 ],'LineWidth',1) % r for red
lgd=legend(NMI_sim_header(2:end));
xlabel('node'); ylabel('LFR CC - group average binary CC');
title('Accuracy of CC Value after Thresholding (Absolute Value)');
lgd.FontSize = 16;
ax = gca; 
ax.FontSize = 16; 
