function [HA_filter]= HA_Filter_KM(HA, LV_mask, Mask_Depth,display)


    HA_filter=HA;
    
    Local_mask=LV_mask;
     
    ListHA=HA(Local_mask>0);
    
    [row,col] = find(Local_mask>0);
    ListDist=Mask_Depth(Local_mask>0);

    ListDist(isnan(ListHA))=nan;
    ListHA(isnan(ListDist))=nan;
    % n=find(isnan(ListHA));
    
    if display
    figure
    scatter(ListDist,ListHA,[],'r','filled');
    axis([0 1 -1.2*100 1.2*100])
    set(gca,'XTick',[0 0.25 0.5 0.75 1]);
    set(gca,'XTickLabel',{'ENDO','','MID','','EPI'});
    set(gca,'YTick',[-90 -45 0 45 90]);
    
    xlabel('')
    ylabel('HA(°)')
    set(gca, 'box', 'off') % remove top x-axis and right y-axis
    set(gcf, 'color', [1 1 1]);
    set(gca, 'color', [1 1 1]);
    ax = gca;
    ax.XColor = 'black';
    ax.YColor = 'black';
    ax.FontSize=15;
    ax.FontWeight='bold';

    legend('off');
    grid off
    end
   

    %row(n)=[];
    %col(n)=[];
    ListDist2=ListDist;
    ListHA2=ListHA;
    
    ListDist2(isnan(ListDist))=[];
    ListHA2(isnan(ListHA))=[];
    
    %% HA has to be flipped
    if nanmedian(ListHA2(find(ListDist2<0.1)))<0
        ListHA=-ListHA;
    end
    
    f = fittype('a*x+b'); 
    fit1 = fit(ListDist2,ListHA2,f,'StartPoint',[1 1]);
    fdata = feval(fit1,ListDist); 
    I = abs(fdata - ListHA) > 1.5*nanstd(ListHA); 
    outliers = excludedata(ListDist,ListHA,'indices',I);
     
    if display
    figure
    scatter(ListDist(~outliers),ListHA(~outliers),[],'r','filled');
    hold on
    scatter(ListDist(outliers),ListHA(outliers),[],'m','filled');
    plot(ListDist,fdata,'-k','LineWidth',4)
    plot((0:0.01:1),fit1(0:0.01:1)+1.5*nanstd(ListHA),':k','LineWidth',4)
    plot((0:0.01:1),fit1(0:0.01:1)-1.5*nanstd(ListHA),':k','LineWidth',4)
    axis([0 1 -1.2*100 1.2*100])
    set(gca,'XTick',[0 0.25 0.5 0.75 1]);
    set(gca,'XTickLabel',{'ENDO','','MID','','EPI'});
    set(gca,'YTick',[-90 -45 0 45 90]);
    
    xlabel('')
    ylabel('HA(°)')
    set(gca, 'box', 'off') % remove top x-axis and right y-axis
    set(gcf, 'color', [1 1 1]);
    set(gca, 'color', [1 1 1]);
    ax = gca;
    ax.XColor = 'black';
    ax.YColor = 'black';
    ax.FontSize=15;
    ax.FontWeight='bold';


    legend('off');
    grid off
    end
    ListOutliers=ListHA(outliers);
    ListDistOutliers=ListDist(outliers);
    ListOutliers(ListDistOutliers<0.5&ListOutliers<0)=-ListOutliers(ListDistOutliers<0.5&ListOutliers<0);
    ListOutliers(ListDistOutliers>0.5&ListOutliers>0)=-ListOutliers(ListDistOutliers>0.5&ListOutliers>0);
    ListHA(outliers)=ListOutliers;

    if display
    figure
    scatter(ListDist(~outliers),ListHA(~outliers),[],'r','filled');
    hold on
    scatter(ListDist(outliers),ListHA(outliers),[],'m','filled');
    axis([0 1 -1.2*100 1.2*100])
    set(gca,'XTick',[0 0.25 0.5 0.75 1]);
    set(gca,'XTickLabel',{'ENDO','','MID','','EPI'});
    set(gca,'YTick',[-90 -45 0 45 90]);
    
    xlabel('')
    ylabel('HA(°)')
    set(gca, 'box', 'off') % remove top x-axis and right y-axis
    set(gcf, 'color', [1 1 1]);
    set(gca, 'color', [1 1 1]);
    ax = gca;
    ax.XColor = 'black';
    ax.YColor = 'black';
    ax.FontSize=15;
    ax.FontWeight='bold';

    legend('off');
    grid off
    end 
    HA_filter(sub2ind(  size(HA_filter), row,col))=ListHA;
    
end
