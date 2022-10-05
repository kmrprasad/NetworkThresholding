%% No Negatives

load('Results_noNegatives.mat');

y1= d_LFRb - d_FCsim_abs;
y2 = d_LFRb - d_FC_perc;
y3 = d_LFRb - d_FC_objL;
y4 = d_LFRb - d_FC_stat;
y5 = d_LFRb - d_FC_tree; 

figure; hold on;
plot(y1,'Color',[1 .8 .8],'LineWidth',5) % r for red
plot(y2,'Color',[.6, .6, .6],'LineWidth',4) % r for red
plot(y3,'Color',[1 .4 .4],'LineWidth',3) % r for red
plot(y4,'Color',[.2 .2 .2],'LineWidth',2) % r for red
plot(y5,'Color',[ .4 .4 1 ],'LineWidth',1) % r for red
lgd=legend(NMI_sim_header(2:end));
xlabel('node'); ylabel('LFR CC - group average binary CC');
title('Accuracy of Density Value after Thresholding (No Negatives)');
lgd.FontSize = 16;
ax = gca; 
ax.FontSize = 16; 


%% Absolute Value

load('Results_absValue.mat');

y1= d_LFRb - d_FCsim_abs;
y2 = d_LFRb - d_FC_perc;
y3 = d_LFRb - d_FC_objL;
y4 = d_LFRb - d_FC_stat;
y5 = d_LFRb - d_FC_tree; 

figure; hold on;
plot(y1,'Color',[1 .8 .8],'LineWidth',5) % r for red
plot(y2,'Color',[.6, .6, .6],'LineWidth',4) % r for red
plot(y3,'Color',[1 .4 .4],'LineWidth',3) % r for red
plot(y4,'Color',[.2 .2 .2],'LineWidth',2) % r for red
plot(y5,'Color',[ .4 .4 1 ],'LineWidth',1) % r for red
lgd=legend(NMI_sim_header(2:end));
xlabel('node'); ylabel('LFR CC - group average binary CC');
title('Accuracy of Density Value after Thresholding (Absolute Value)');
lgd.FontSize = 16;
ax = gca; 
ax.FontSize = 16; 
