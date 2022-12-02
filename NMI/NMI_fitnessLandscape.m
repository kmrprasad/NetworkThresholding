function NMI_fitnessLandscape(DCN,FCN,optima_threshold_d,optima_threshold_f,thresholds,M0_d,M1_d,M0_f,M1_f)
%NMI_FITNESSLANDSCAPE Summary of this function goes here
%   Detailed explanation goes here

header={'Density & Degree'; 'Density & Degree';...
    'Lambda'; 'Transitivity';
    'Average CC'; 'Efficiency';...
    'Modularity';'Average BC';'Assortativity'};

numThresh=length(thresholds);
numSubs=size(DCN,3);

F=zeros(numThresh,numThresh,numSubs);

X=thresholds;
Y=thresholds;

for n=1:numSubs
    for x=1:numThresh
        for y=1:numThresh
                DCON_bin=ThresholdNetwork(DCN(:,:,n), thresholds(x));
                FCON_bin=ThresholdNetwork(FCN(:,:,n), thresholds(y));

                try
                    [communityNMI, ~, ~, ~, ~] = NMI_NetMod_Main(DCON_bin,FCON_bin);
                catch
                   communityNMI=0; 
                end
                
                F(x,y,n)=communityNMI;
        end
    end
end

Z=mean(F,3);

Z(isnan(Z))=0;

figure;
[M,c]=contour(X,Y,Z,20);
c.LineWidth = 3;
colormap('hot');
colorbar;
hold on;
for g = 1:size(optima_threshold_d,1)
    xi=mean(optima_threshold_f(g,:));
    yi=mean(optima_threshold_d(g,:));
    scatter(xi,yi,'m','filled');
    c=header{g};
    text(xi, yi, c, 'Fontsize', 16);
end

yMin=mean(thresholds(M0_d));
yMax=mean(thresholds(M1_d));
xMin=mean(thresholds(M0_f));
xMax=mean(thresholds(M1_f));

yline(yMin);
yline(yMax);

xline(xMin);
xline(xMax);

title('Fitness Landscape');
ylabel('DCN threshold');
xlabel('FCN threshold');

xBuffer=0.1*(xMax-xMin);
yBuffer=0.1*(yMax-yMin);

xlim([xMin-xBuffer, xMax+xBuffer]);
ylim([yMin-yBuffer,yMax+yBuffer]);


end





