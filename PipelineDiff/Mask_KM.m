function [Dcm2]= Mask_KM(Dcm,min, max)

% Generate a mask matrix  
%
% SYNTAX:  [Dcm2]= Mask_KM(Dcm,min, max)
%  
%
% INPUTS:   Dcm - Image matrix
%                 [y x slices ..]
%
%           min - Seuil minimal for the mask generation
%           
%           max - Seuil maximal for the mask generation 
%
%        
% OUTPUTS:  Dcm2 - image matrix 
%                 [y x slices ..]
%
%
% Kevin Moulin 01.13.2020
% Kevin.Moulin.26@gmail.com
% Ennis Lab @ UCLA: http://mrrl.ucla.edu
% Ennis Lab @ Stanford: https://med.stanford.edu/cmrgroup/software.html

    Dcm2=[];
    disp('Mask creation') 
    Dcm2=Dcm;
    
    Dcm2(Dcm2>max)=0;
    Dcm2(Dcm2<min)=0;

end