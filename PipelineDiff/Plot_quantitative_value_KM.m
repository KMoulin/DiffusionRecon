
function Plot_quantitative_value_KM(TD,Value,Color,RangeX,RangeY)

%% rHA_range vs TD all 
hold on

ln = plot([0.5 0.5], [-150 150],'-'); hold on;
ln.Color = [.5 .5 .5];
ln.LineWidth = 1;                

ln = plot([0 1], [0 0],'-'); hold on;
ln.Color = [.5 .5 .5];
ln.LineWidth = 1;            


Bplot = plot(TD',Value', '--','LineWidth',2,'Color',[0.5 0.5 0.5]);
for cpt=1:1:length(TD)
    
    Bplot = scatter(TD(cpt),Value(cpt),100,'o','MarkerFaceAlpha',7/8);
   % Bplot.MarkerSize = 20;
    Bplot.MarkerFaceColor = Color (cpt,:);
    Bplot.MarkerEdgeColor = 0.5*Color(cpt,:);
   
    
end
axis([RangeX(1) RangeX(end) RangeY(1) RangeY(end)])
set(gca,'XTick',RangeX);
xlabel('')
ylabel('')
set(gca, 'box', 'on') % remove top x-axis and right y-axis
set(gcf, 'color', [1 1 1]);
set(gca, 'color', [1 1 1]);
ax = gca;
ax.XColor = 'black';
ax.YColor = 'black';
ax.FontSize=25;
ax.FontWeight='bold';

ax.LineWidth = 1.5;
ax.Box=0;
ax.YTick=RangeY;

legend('off');
grid off

end