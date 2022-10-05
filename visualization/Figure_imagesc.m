%% No Negatives

load('Results_absValue.mat');
figure;
subplot(1,2,1)
imagesc(LFR_orig);
colormap('pink');
colorbar;
title('LFR network');
subplot(1,2,2);
imagesc(mean(FC_sim,3));
title('Average FCsim');
colormap('pink');
colorbar;