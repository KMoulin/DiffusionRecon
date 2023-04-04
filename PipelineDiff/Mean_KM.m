function [Dcm2 enum2]= Mean_KM(Dcm, enum)

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
% Kevin Moulin 08.14.2017
% Kevin.Moulin.26@gmail.com
% Ennis Lab @ UCLA; http://mrrl.ucla.edu

    enum2=enum; 
    Dcm2=[];
    disp('Mean data') 
    h = waitbar(0,'Mean data...');
    
    hf = 1/3*ones(3,1);
    H = hf*hf';


    for cpt_set=1:1:enum.nset
        for cpt_slc=1:1:enum.datasize(cpt_set).slc
         for cpt_b=1:1:enum.datasize(cpt_set).b     
           for cpt_dir=1:1: enum.dataset(cpt_set).slc(cpt_slc).b(cpt_b).nb_dir     
               for cpt_avg=1:1: enum.dataset(cpt_set).slc(cpt_slc).b(cpt_b).dir(cpt_dir).nb_avg
                    tmp = filter2(H,squeeze((Dcm(:,:,cpt_slc,cpt_b,cpt_dir,cpt_avg,cpt_set))));
                    Dcm2(:,:,cpt_slc,cpt_b,cpt_dir,cpt_avg,cpt_set)=tmp;                   
               end
           end
         end   
         waitbar(cpt_slc/enum.datasize(cpt_set).slc,h);
        end
    end
    close(h);    

end