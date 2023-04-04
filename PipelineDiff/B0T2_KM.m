function [T2Map T0Map]= B0T2_KM(DcmB0, enum)


%  Generate T2 maps: 
%  
%
% SYNTAX:  [ADC]= ADC_KM(Dcm, DcmB0, enum)
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
%          
% OUTPUTS:  ADC - ADC image matrix
%                 [y x slices b-values directions averages]
%
%
% Kevin Moulin 08.14.2017
% Kevin.Moulin.26@gmail.com
% Ennis Lab @ UCLA; http://mrrl.ucla.edu

    T2Map=[];
    T0Map=[];
    disp('T2 calculation') 
   % h = waitbar(0,'T2 calculation...');
    for cpt_slc=1:1:size(enum.slc,2)              
        [tmpT2Map tmpT0Map]= T2Fit_KM(squeeze(DcmB0(:,:,cpt_slc,:)), enum.TE);  
        T2Map(:,:,cpt_slc)=tmpT2Map;
        T0Map(:,:,cpt_slc)=tmpT0Map;
        %waitbar(cpt_slc/size(enum.slc,2),h);
    end
    %close(h);    
end