function [Dcm2 enum2]= Average_KM(Dcm, enum)

% Perform averaging on the DWI and nDWI matrices
% 
% Update the number of averages in enum
%
%
% SYNTAX:  [Dcm2 enum2]= Average_KM(Dcm, enum);
%  
%
% INPUTS:   Dcm - DWI image matrix
%                 [y x slices b-values directions averages dataset]
%           
%           enum - Structure which contains information about the dataset 
%          
% OUTPUTS:  Dcm2 - DWI image matrix 
%                 [y x slices b-values directions averages dataset]
%           
%           enum2 - Structure which contains information about the dataset 
%           
%
% Kevin Moulin 01.13.2020
% Kevin.Moulin.26@gmail.com
% Ennis Lab @ UCLA: http://mrrl.ucla.edu
% Ennis Lab @ Stanford: https://med.stanford.edu/cmrgroup/software.html

    enum2=enum; 
    Dcm2=[];
    disp('Average data') 
    h = waitbar(0,'Average data...');
    for cpt_set=1:1:enum.nset
        for cpt_slc=1:1:enum.datasize(cpt_set).slc
         for cpt_b=1:1:enum.datasize(cpt_set).b     
           for cpt_dir=1:1: enum.dataset(cpt_set).slc(cpt_slc).b(cpt_b).nb_dir                                 
               Dcm2(:,:,cpt_slc,cpt_b,cpt_dir,1,cpt_set)=mean(Dcm(:,:,cpt_slc,cpt_b,cpt_dir,1:enum.dataset(cpt_set).slc(cpt_slc).b(cpt_b).dir(cpt_dir).nb_avg,cpt_set),6);              
               enum2.dataset(cpt_set).slc(cpt_slc).b(cpt_b).dir(cpt_dir).nb_avg=1;
           end
         end   
         waitbar(cpt_slc/enum.datasize(cpt_set).slc,h);
        end
    end
    close(h);    

end