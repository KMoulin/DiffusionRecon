function [Mask_AHA] = ROI2AHA_KM (Dcm, P_Endo, P_Epi)

% Generate a mask matrix with 6 segments corresponding to the AHA cardiac segmentation  
%
% SYNTAX:  [Mask_AHA] = ROI2AHA_KM (Mask, P_Endo, P_Epi)
%  
% INPUTS:   Dcm - Image matrix
%                 [y x slices]
%
%           P_Endo - List of Coordinates of the Endocardium ROI
%                 
%           P_Epi - List of Coordinates of the Endocardium ROI
%           
% OUTPUTS:  Mask_AHA - Mask matrix 
%                 [y x slices AHA_segments]
%
%
% Kevin Moulin 08.14.2017
% Kevin.Moulin.26@gmail.com
% Ennis Lab @ UCLA; http://mrrl.ucla.edu

        for cpt_slc=1:1:size(Dcm,3)
            figure (99)
            imagesc(Dcm(:,:,cpt_slc,1,1,1,1,1,1)); 
            title('Choose the jonction between RV/LV')       
            [px py] = ginput(1);

            center=[ mean( P_Endo(:,1,cpt_slc) );mean( P_Endo(:,2,cpt_slc) )];
            center2=[ mean( P_Epi(:,1,cpt_slc) );mean( P_Epi(:,2,cpt_slc) )];
            Mask_AHA(:,:,cpt_slc,:) = Divide_n_Rule(center,[px py], squeeze(Dcm(:,:,cpt_slc,1,1,1,1,1,1)));
            close(99)
        end
        
        Mask_AHA(Mask_AHA>0)=1;
end


function [Mask_AHA] = Divide_n_Rule(Center,VDVG, Mask)


Mask_AHA = zeros(size(Mask,1), size(Mask,2) , 6);



VectVDVG(1)=(VDVG(1)-Center(1));
VectVDVG(2)=(VDVG(2)-Center(2));
VectVDVG(3)=0;
Mask_Tmp=[];
for x = 1:1:size(Mask,2)
    for y = 1:1:size(Mask,1)
        VectPix(1)=(x-Center(1));
        VectPix(2)=(y-Center(2));
        VectPix(3)=0;

       % Unsigned angle between the two vectors
        theta = acos(dot(VectVDVG / norm(VectVDVG), VectPix / norm(VectPix)));

        % Determine the sign of the angle
        sgn = sign(cross(VectVDVG, VectPix));

        % Apply the sign and use mod to make it between 0 and 2*pi
         Mask_Tmp(y,x) = mod(theta * (-1)^(sgn(3) < 0), 2*pi);
         
      if  Mask_Tmp(y,x)>0 && 1*2*pi/6 >= Mask_Tmp(y,x)
        Mask_AHA(y,x,1) = Mask(y,x);
      elseif Mask_Tmp(y,x)>1*2*pi/6 && 2*2*pi/6 >= Mask_Tmp(y,x)
        Mask_AHA(y,x,2) = Mask(y,x);
      elseif Mask_Tmp(y,x)>2*2*pi/6 && 3*2*pi/6 >= Mask_Tmp(y,x)
        Mask_AHA(y,x,3) = Mask(y,x);  
      elseif Mask_Tmp(y,x)>3*2*pi/6 && 4*2*pi/6 >= Mask_Tmp(y,x)
        Mask_AHA(y,x,4) = Mask(y,x);  
      elseif Mask_Tmp(y,x)>4*2*pi/6 && 5*2*pi/6 >= Mask_Tmp(y,x)
        Mask_AHA(y,x,5) = Mask(y,x);  
      else
        Mask_AHA(y,x,6) = Mask(y,x);  
      end
          
    end
end


end
% 