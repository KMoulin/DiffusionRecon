function [RA]= Plot_Wall_Depth_KM( HA,Mask_Depth)

    Local_mask=Mask_Depth;
    ListHA=HA(~isnan(Local_mask));
    ListDist=Mask_Depth(~isnan(Local_mask));
    WD_pitch= 1/20;
    WD_all = WD_pitch:WD_pitch:1; %wall depth bins
    colors=[1,0,0]; 
    figure,
    scatter(ListDist,ListHA,[],colors,'filled'); hold on;

    % Bin and plot HA confidence intervals
     cpt = 1;
    for WD = WD_all  

        I = find(ListDist < WD & ListDist >= WD - WD_pitch);

        tmpDat = ListHA(I);
        Med(cpt) = median(tmpDat);
        tmpRange = quantile(tmpDat,[0.25 0.75]);

        x = WD - 0.025;


        plot(x,Med(cpt),'k+','MarkerSize',20); hold on; %plot bin medians
        plot([x-0.01 x+0.01], [tmpRange(1) tmpRange(1)], 'k-','MarkerSize',20,'LineWidth',3);        %plot lower CI
        plot([x-0.01 x+0.01], [tmpRange(2) tmpRange(2)], 'k-','MarkerSize',20,'LineWidth',3);        %plot upper CI
        cpt = cpt + 1;

    end

    RA=Med(end)-Med(1);
    
    axis([0 1 -100 100]);



end