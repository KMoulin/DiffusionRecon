function Overlay_img_KM(im1,im2,mask,color1,color2,range1,range2,type)



if type==[]
   type=1; 
end


for cpt_mask=1:1:size(mask,3)
    Alpha=zeros(size(mask));
    im1(mask~=0)=0;
    im2(mask==0)=0;
    Alpha(mask==0)=1;

    figure,
    ax1 = axes;
    imagesc(im1(:,:,cpt_mask),'AlphaData',Alpha(:,:,cpt_mask))
    ax1.CLim=[range1(1) range1(2)];
    ax1.Visible = 'off';
    ax1.XTick = [];
    ax1.YTick = [];
    if type==1
        axis image;
    elseif type==2
        axis square;
    elseif type==3
        axis equal;
    end
    view(2)
    ax2 = axes;
    if  size(im2,3)==3
         imagesc(squeeze(im2(:,:,:)),'AlphaData',1-Alpha(:,:,cpt_mask))
    elseif size(im2,4)==1
        imagesc(im2(:,:,cpt_mask),'AlphaData',1-Alpha(:,:,cpt_mask))
    else
        imagesc(squeeze(im2(:,:,cpt_mask,:)),'AlphaData',1-Alpha(:,:,cpt_mask))
    end
    ax2.CLim=[range2(1) range2(2)];
    ax2.Visible = 'off';
    ax2.XTick = [];
    ax2.YTick = [];
    if type==1
        axis image;
    elseif type==2
        axis square;
    elseif type==3
        axis equal;
    end

    linkaxes([ax1,ax2])

    colormap(ax1,color1)
    colormap(ax2,color2)    
end
end