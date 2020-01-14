function [Dcm2]= Norm_KM(Dcm, enum)
    
%  Normalize DWI matrices by dividing it by the nDWI matrices 
%
%
% SYNTAX:  [Dcm2]= Norm_KM(Dcm, enum)
%  
%
% INPUTS:   Dcm - DWI image matrix
%                 [y x slices b-values directions averages dataset]
%           
%           enum - Structure which contains information about the dataset 
%
%          
% OUTPUTS:  Dcm2 - DWI image matrix
%                 [y x slices b-values directions averages dataset]
%
%
% Kevin Moulin 01.13.2020
% Kevin.Moulin.26@gmail.com
% Ennis Lab @ UCLA: http://mrrl.ucla.edu
% Ennis Lab @ Stanford: https://med.stanford.edu/cmrgroup/software.html

    Dcm2=[];
    disp('Norm data (S/S0)') 
    h = waitbar(0,'Norm (S/S0)...');
    for cpt_set=1:1:enum.nset
        for cpt_slc=1:1:enum.datasize(cpt_set).slc
            for cpt_b=1:1:enum.datasize(cpt_set).b     
              for cpt_dir=1:1: enum.dataset(cpt_set).slc(cpt_slc).b(cpt_b).nb_dir           
                 for cpt_avg=1:1: enum.dataset(cpt_set).slc(cpt_slc).b(cpt_b).dir(cpt_dir).nb_avg                  
                     Dcm2(:,:,cpt_slc,cpt_b,cpt_dir,cpt_avg,cpt_set)=squeeze(Dcm(:,:,cpt_slc,cpt_b,cpt_dir,cpt_avg,cpt_set))./squeeze(Dcm(:,:,cpt_slc,1,1,1,cpt_set));                     
                 end         
                end
            end
           waitbar(cpt_slc/enum.datasize(cpt_set).slc,h);
        end
    end
    close(h);    
    Dcm2(isnan(Dcm2))=0;
    Dcm2(isinf(Dcm2))=0;
end