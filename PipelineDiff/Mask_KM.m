function [Dcm2]= Mask_KM(Dcm,min, max)

% Generate a mask matrix  
%
% SYNTAX:  [Dcm2]= Mask_KM(Dcm,min, max)
%  
%
% INPUTS:   Dcm - Image matrix
%                 [y x slices b-values directions averages dataset]
%
%           min - Seuil minimal for the mask generation
%           
%           max - Seuil maximal for the mask generation 
%
%        
% OUTPUTS:  Dcm2 - image matrix 
%                 [y x slices b-values directions averages dataset]
%
%
% Kevin Moulin 08.14.2017
% Kevin.Moulin.26@gmail.com
% Ennis Lab @ UCLA; http://mrrl.ucla.edu

    Dcm2=[];
    disp('Mask creation') 
    Dcm2=Dcm;
    
    Dcm2(Dcm2>max)=0;
    Dcm2(Dcm2<min)=0;

end