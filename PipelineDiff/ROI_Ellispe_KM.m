function[Dcm2,DcmROI,Mask,MaskROI]= ROI_Ellispe_KM(Dcm,nROI)

%  Draw a ROI for each Slice and apply it to the DWI matrix 
%  
% SYNTAX:  [Dcm2,Dcm3,Mask,Mask_mult]= ROI_Ellispe_KM(Dcm,nROI)
%  
%
% INPUTS:   Dcm - Image matrix
%                 [y x slices]
%
%           nROI - Number of ROI wanted
%                 [y x slices]
%        
% OUTPUTS:  Dcm2 - Image matrix after application of all the ROIs
%                 [y x slices]
%
%           DcmROI - Image matrix for each ROI 
%                 [y x slices .. nROI]
%
%           Mask - Mask matrix 
%                 [y x slices]
%
%           MaskROI - Mask matrix for each ROI 
%                 [y x slices .. nROI]
%
%
% Kevin Moulin 01.13.2020
% Kevin.Moulin.26@gmail.com
% Ennis Lab @ UCLA: http://mrrl.ucla.edu
% Ennis Lab @ Stanford: https://med.stanford.edu/cmrgroup/software.html

    Dcm2=[];
    DcmROI=[];
    MaskROI=[];
    Mask=[];
    
    disp('Create ROI') 
    for cpt_tube=1:1:nROI
       position=[];
       h = imagesc(Dcm(:,:,1,1,1,1));
       h = imellipse;
       position = wait(h);
       for cpt_slc=1:1:size(Dcm,3)
        MaskROI(:,:,cpt_slc,cpt_tube) = poly2mask(position(:,1),position(:,2),size(Dcm,1),size(Dcm,2));
       end
    end
    
    Mask=max(MaskROI,[],4);
    MaskROI(MaskROI>0)=1;
    Mask(Mask>0)=1;
    for cpt_slc=1:1:size(Dcm,3)
     for cpt_b=1:1:size(Dcm,4)
      for cpt_dir=1:1:size(Dcm,5)
       for cpt_avg=1:1:size(Dcm,6)  
          Dcm2(:,:,cpt_slc,(cpt_b),cpt_dir,cpt_avg)=squeeze(Dcm(:,:,cpt_slc,(cpt_b),cpt_dir,cpt_avg)).*Mask(:,:,cpt_slc);
          for cpt_tube=1:1:nROI
                DcmROI(:,:,cpt_slc,(cpt_b),cpt_dir,cpt_avg,cpt_tube)=squeeze(Dcm(:,:,cpt_slc,(cpt_b),cpt_dir,cpt_avg)).*MaskROI(:,:,cpt_slc,cpt_tube);
          end
       end
      end
     end
    end
end