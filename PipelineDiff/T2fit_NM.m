function [T2 S0]= T2fit_NM(S1, S2, enum)

% Perform averaging on the DWI and nDWI matrices
% 
% Update the number of averages in enum
%
%
% SYNTAX:  [Dcm2 DcmB02 enum2]= Average_KM(Dcm, DcmB0, enum);
%  
%
% INPUTS:   Dcm - DWI image matrix
%                 [y x slices b-values directions averages]
%
%           DcmB0 - Non diffusion weighted image matrix 
%                  [y x slices b-values directions averages]
%           
%           enum - Structure which contains information about the dataset 
%          
% OUTPUTS:  Dcm2 - DWI image matrix 
%                 [y x slices b-values directions averages]
%
%           DcmB02 - Non diffusion weighted image matrix 
%                  [y x slices b-values directions averages]
%           
%           enum2 - Structure which contains information about the dataset 
%
% Kevin Moulin 08.14.2017
% Kevin.Moulin.26@gmail.com
% Ennis Lab @ UCLA; http://mrrl.ucla.edu


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