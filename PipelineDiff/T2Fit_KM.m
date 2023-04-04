function [T2 S0]= T2fit_KM(S1, S2, enum)

% Generate T2 maps and the S0 maps from Two differents TE images
%
% SYNTAX:  [T2 S0]= T2fit_KM(S1, S2, enum);
%  
%
% INPUTS:   S1 - Long TE matrix images corresponding to enum.TE(1)
%                 [y x slices]
%
%           S2 - Short TE matrix images corresponding to enum.TE(2) 
%                  [y x slices]
%           
%           enum - Structure which contains information about the dataset 
%          
% OUTPUTS:  T2 - T2 image matrix 
%                 [y x slices]
%
%           S0 - S0 matrix
%                 [y x slices [S01 S02]]           
%
%
% Kevin Moulin 01.13.2020
% Kevin.Moulin.26@gmail.com
% Ennis Lab @ UCLA: http://mrrl.ucla.edu
% Ennis Lab @ Stanford: https://med.stanford.edu/cmrgroup/software.html


%
%
% S2=exp(-TE2/T2) 
% 
% S1=exp(-TE1/T2)
% 
% S2/S1=exp(-TE2/T2 + TE1/2)
% 
% log(S2/S1) = (TE1-TE2)/T2 
% 
% T2=(TE1-TE2)/log(S2/S1)
    
T2=[];
    disp('T2 fit') 
    h = waitbar(0,'T2 fit...');
    for cpt_slc=1:1:size(enum.slc,2)                 
       T2(:,:,cpt_slc)=(enum.TE(1)-enum.TE(2))./log(squeeze(S2(:,:,cpt_slc,1,1,1))./squeeze(S1(:,:,cpt_slc,1,1,1)));
       S0(:,:,cpt_slc,1)=squeeze(S1(:,:,cpt_slc,1,1,1))./squeeze(exp(-enum.TE(1)./T2(:,:,cpt_slc)));
       S0(:,:,cpt_slc,2)=squeeze(S2(:,:,cpt_slc,1,1,1))./squeeze(exp(-enum.TE(2)./T2(:,:,cpt_slc)));
       waitbar(cpt_slc/size(enum.slc,2),h);
    end
    
    T2(isinf(T2))=0;
    T2(isnan(T2))=0;
    close(h);    
end