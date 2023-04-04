function[Dcm2,Dcm3,mask,mask_mult]= ROI_Ellispe_Tube_KM(Dcm,nTube)

    Dcm2=[];
    Dcm3=[];
    mask_mult=[];
    mask=[];
    
    disp('Create ROI') 
    for cpt_tube=1:1:nTube
       position=[];
       h = imagesc(Dcm(:,:,1,1,1,1));
       h = imellipse;
       position = wait(h);
       mask_mult(:,:,cpt_tube) = poly2mask(position(:,1),position(:,2),size(Dcm,1),size(Dcm,2));
    end
    mask=max(mask_mult,[],3);
      

      Dcm2(:,:)=squeeze(Dcm(:,:).*mask);
      for cpt_tube=1:1:nTube
            Dcm3(:,:,cpt_tube)=squeeze(Dcm(:,:).*mask_mult(:,:,cpt_tube));
      end

end