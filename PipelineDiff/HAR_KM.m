function [HAR, Med]= HAR_KM(ListHA, ListDist,Range,Color)


if nargin==2
Range=[-90 90];
Color=[1 0 0];
bDisplay=0;
elseif nargin==3
Color=[1 0 0];
bDisplay=1;
elseif nargin==4
    bDisplay=1;
else
 error('Need at least ListHA and ListDist');   
end
    
ImgParam.edges=[Range(1):Range(2)/2:Range(2)]; 
ImgParam.BoxCenters=[2, 6, 10, 14, 18];
                                                                                                                                
ImgParam.py=[40:200];
ImgParam.px=[30:200];

ImgParam.nbins = 9;
ImgParam.WD_pitch= 1/ImgParam.nbins;
ImgParam.WD_all = ImgParam.WD_pitch:ImgParam.WD_pitch:1; %wall depth bins

Med = [];

Med_all = [];
Med_all = zeros(3,ImgParam.nbins)*nan;

ImgParam.colors(1,:)=[0;158/255;115/255]; 
ImgParam.colors(2,:)=[230/255;159/255;0/255]; 
ImgParam.colors(3,:)=[86/255;180/255;233/255]; 

ImgParam.colors(1,:)=Color;

scanNum=1;

%figure('Position',[500 500 400 300])
if bDisplay
ln = plot([0.5 0.5], [-150 150],'-'); hold on;
ln.Color = [.5 .5 .5];
ln.LineWidth = 1;                

ln = plot([0 1], [0 0],'-'); hold on;
ln.Color = [.5 .5 .5];
ln.LineWidth = 1;

scatter(ListDist,ListHA,[],ImgParam.colors(scanNum,:),'filled');
end
% Bin and plot HA confidence intervals
 count = 1;
for WD = ImgParam.WD_all  

    I = find(ListDist < WD & ListDist >= WD - ImgParam.WD_pitch);

    tmpDat = ListHA(I);
    Med(count) = nanmedian(tmpDat);
    tmpRange = quantile(tmpDat,[0.25 0.75]);

    x = WD - ImgParam.WD_pitch/2;
    if bDisplay
        crossplot = plot(x,Med(count),'k+','MarkerSize',20); hold on; %plot bin medians
        crossplot.LineWidth = 2;
        plot([x-0.02 x+0.02], [tmpRange(1) tmpRange(1)], 'k-','MarkerSize',20,'LineWidth',4);        %plot lower CI
        plot([x-0.02 x+0.02], [tmpRange(2) tmpRange(2)], 'k-','MarkerSize',20,'LineWidth',4);        %plot upper CI
    end
    count = count + 1;

end
if bDisplay
    axis([0 1 1.2*min(ImgParam.edges) 1.2*max(ImgParam.edges)])
    set(gca,'XTick',[0 0.25 0.5 0.75 1]);
    set(gca,'XTickLabel',{'ENDO','','MID','','EPI'});
    set(gca,'XTick',[]);

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

    ax.YTick=ImgParam.edges;

    legend('off');
    grid off;
    hold off;
end
HAR=abs(Med(end)-Med(1));

end