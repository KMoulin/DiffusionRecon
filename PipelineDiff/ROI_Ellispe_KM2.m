function[Mask_mult]= ROI_Ellispe_KM2(Dcm,nROI)

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
       
        Mask_mult(:,:,cpt_tube) = poly2mask(position(:,1),position(:,2),size(Dcm,1),size(Dcm,2));
    end
end