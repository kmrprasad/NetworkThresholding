function [Tree] = MaxTreeFig(Adj,coordinates)
%ControlCentor.... from %https://www.mathworks.com/matlabcentral/fileexchange/23276-maximum-weight-spanning-tree-undirected
% Adj can be considered as the cost matrix
% Tree is the matrix saving the tree, and Cost is summation of all cost
% upon arcs.
[ Tree,~ ] =  UndirectedMaximumSpanningTree ( Adj );
%h1 = view(biograph( Tree ));
figure;
imagesc(Tree);
% 3d plot
Stree=sparse(Tree);
figure;
%
subplot(1,3,1);
nYZ=[coordinates{:,3};coordinates{:,4}]';
nYZ(:,1)=-nYZ(:,1);
gplot(Stree,nYZ,'-b');
hold on;
scatter(nYZ(:,1),nYZ(:,2),'o','red','filled'); set(gca,'XTickLabel',[]); set(gca,'YTickLabel',[]);
%
subplot(1,3,2);
XY=[coordinates{:,2};coordinates{:,3}]';
gplot(Stree,XY,'-b');
hold on;
scatter(XY(:,1),XY(:,2),'o','red','filled'); set(gca,'XTickLabel',[]); set(gca,'YTickLabel',[]);
% XZ=[coordinates{:,2};coordinates{:,4}]';
% gplot(Stree,XZ);
% hold on;
% scatter(XZ(:,1),XZ(:,2),'o','red','filled');
%
subplot(1,3,3);
YZ=[coordinates{:,3};coordinates{:,4}]';
gplot(Stree,YZ,'-b');
hold on;
scatter(YZ(:,1),YZ(:,2),'o','red','filled'); set(gca,'XTickLabel',[]); set(gca,'YTickLabel',[]);

figure;
Gt=graph(Tree);
plot(Gt,'XData',cell2mat(coordinates(:,2)),'YData',cell2mat(coordinates(:,3)),'ZData',cell2mat(coordinates(:,4)));
%scatter3(coordinates{:,2},coordinates{:,3},coordinates{:,4},'o','red','filled')
end

