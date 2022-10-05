function  BarGraphs(Data,GraphTitle,GroupNames)

data_mean=mean(Data,1);
data_std=std(Data,[],1);


figure;
y=data_mean;
errorhigh=data_std;
x=1:length(y);
b=bar(x,y);               
b.FaceColor=[0.7,0.7,0.7];
hold on


title(GraphTitle);
xticks(x);
set(gca, 'XTickLabel',GroupNames,'FontSize', 16);
xtickangle(45);
ylabel('NMI');
er = errorbar(x,y,errorhigh,'vertical');    
er.Color = [0 0 0];            
er.LineStyle = 'none';  

ty=y+(1.3*errorhigh);
formatSpec = '%.3f';

z=string(y);
for baby_y=1:length(y)
   z(baby_y)=num2str(y(baby_y),formatSpec);
end

text(x,ty,z,'HorizontalAlignment','center');

end

