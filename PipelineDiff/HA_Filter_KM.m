function [HA_filter]= HA_Filter_KM(HA, LV_mask, Mask_Depth)


    HA_filter=HA;
    
    Local_mask=LV_mask;
    
    ListHA=HA(Local_mask>0);
    
    [row,col] = find(Local_mask>0);
    ListDist=Mask_Depth(Local_mask>0);

    ListDist(isnan(ListHA))=nan;
    ListHA(isnan(ListDist))=nan;
    n=find(isnan(ListHA));
    ListDist(isnan(ListDist))=[];
    ListHA(isnan(ListHA))=[];

    row(n)=[];
    col(n)=[];
    
    f = fittype('a*x+b'); 
    fit1 = fit(ListDist,ListHA,f,'StartPoint',[1 1]);
    fdata = feval(fit1,ListDist); 
    I = abs(fdata - ListHA) > 1.5*std(ListHA); 
    outliers = excludedata(ListDist,ListHA,'indices',I);
 
    figure
    scatter(ListDist(~outliers),ListHA(~outliers),[],'r','filled');
    hold on
    scatter(ListDist(outliers),ListHA(outliers),[],'m','filled');
     
    ListOutliers=ListHA(outliers);
    ListDistOutliers=ListDist(outliers);
    ListOutliers(ListDistOutliers<0.5&ListOutliers<0)=-ListOutliers(ListDistOutliers<0.5&ListOutliers<0);
    ListOutliers(ListDistOutliers>0.5&ListOutliers>0)=-ListOutliers(ListDistOutliers>0.5&ListOutliers>0);
    ListHA(outliers)=ListOutliers;

    figure
    scatter(ListDist(~outliers),ListHA(~outliers),[],'r','filled');
    hold on
    scatter(ListDist(outliers),ListHA(outliers),[],'m','filled');
    
    HA_filter(sub2ind(  size(HA_filter), row,col))=ListHA;
    
end
