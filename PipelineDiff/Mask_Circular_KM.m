function [Dcm2]= Mask_Circular_KM(Dcm,rad)

% Generate a mask matrix  
%
% SYNTAX:  [Dcm2]= Mask_KM(Dcm,min, max)
%  
%
% INPUTS:   Dcm - Image matrix
%                 [y x slices b-values directions averages]
%
%           min - Seuil minimal for the mask generation
%           
%           max - Seuil maximal for the mask generation 
%
%        
% OUTPUTS:  Dcm2 - image matrix 
%                 [y x slices b-values directions averages]
%
%
% Kevin Moulin 08.14.2017
% Kevin.Moulin.26@gmail.com
% Ennis Lab @ UCLA; http://mrrl.ucla.edu

    Dcm2=[];
    disp('Mask creation') 
    
    Mask = createCirclesMask([size(Dcm,1) size(Dcm,2)] ,[round(size(Dcm,2)/2) round(size(Dcm,1)/2)],rad);
    
    Vect=size(Dcm);
    Vect(1)=1;
    Vect(2)=1;
    Dcm2 = Dcm.*repmat(Mask,Vect);
  
end