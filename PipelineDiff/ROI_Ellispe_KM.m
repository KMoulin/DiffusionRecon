function[Dcm2,Dcm3,Mask,Mask_mult]= ROI_Ellispe_KM(Dcm,nROI)

%  Draw n ROI for the current Slice and apply it to the DWI matrix 
%  
% SYNTAX:  [Dcm2,Dcm3,mask,mask_mult]= ROI_Ellispe_KM(Dcm,nTube)
%  
%
% INPUTS:   Dcm - Image matrix
%                 [y x]
%
%        
% OUTPUTS:  P_Endo - List of Coordinates of the Endocardium ROI
%                 
%           P_Epi - List of Coordinates of the Epicardium ROI
%
%           Dcm2 - DWI image matrix
%           [y x slices b-values directions averages dataset]
%
% Kevin Moulin 08.14.2017
% Kevin.Moulin.26@gmail.com
% Ennis Lab @ UCLA; http://mrrl.ucla.edu

    Dcm2=[];
    Dcm3=[];
    Mask_mult=[];
    Mask=[];
    
    disp('Create ROI') 
    for cpt_tube=1:1:nROI
       position=[];
       h = imagesc(Dcm(:,:,1,1,1,1));
       h = imellipse;
       position = wait(h);
       for cpt_slc=1:1:size(Dcm,3)
        Mask_mult(:,:,cpt_slc,cpt_tube) = poly2mask(position(:,1),position(:,2),size(Dcm,1),size(Dcm,2));
       end
    end
    
    Mask=max(Mask_mult,[],4);
    Mask_mult(Mask_mult>0)=1;
    Mask(Mask>0)=1;
    for cpt_slc=1:1:size(Dcm,3)
     for cpt_b=1:1:size(Dcm,4)
      for cpt_dir=1:1:size(Dcm,5)
       for cpt_avg=1:1:size(Dcm,6)  
          Dcm2(:,:,cpt_slc,(cpt_b),cpt_dir,cpt_avg)=squeeze(Dcm(:,:,cpt_slc,(cpt_b),cpt_dir,cpt_avg)).*Mask(:,:,cpt_slc);
          for cpt_tube=1:1:nROI
                Dcm3(:,:,cpt_slc,(cpt_b),cpt_dir,cpt_avg,cpt_tube)=squeeze(Dcm(:,:,cpt_slc,(cpt_b),cpt_dir,cpt_avg)).*Mask_mult(:,:,cpt_slc,cpt_tube);
          end
       end
      end
     end
    end
end