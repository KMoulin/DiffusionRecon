function Dcm2= T2correction_KM(Dcm,T2Map, enum)

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

 
   % Dcm2=squeeze(Dcm);
    
   % Dcm2(:,:,:,1:(end-enum.nRestore))=Dcm2(:,:,:,1:(end-enum.nRestore)).*repmat(exp(-T2Map/enum.TE(1)),1,1,1,size(Dcm2,4)-enum.nRestore);
   % Dcm2(:,:,:,end-enum.nRestore+1:end)=Dcm2(:,:,:,end-enum.nRestore+1:end).*repmat(exp(-T2Map/enum.TE(2)),1,1,1,enum.nRestore);
  
    
    disp('T2 Correction') 
    h = waitbar(0,'T2 Correction...');
    for cpt_set=1:1:enum.nset
         for cpt_slc=1:1:enum.datasize(cpt_set).slc
            for cpt_b=1:1:enum.datasize(cpt_set).b     
              for cpt_dir=1:1: enum.dataset(cpt_set).slc(cpt_slc).b(cpt_b).nb_dir           
                 for cpt_avg=1:1: enum.dataset(cpt_set).slc(cpt_slc).b(cpt_b).dir(cpt_dir).nb_avg                   
                            Dcm2(:,:,cpt_slc,cpt_b,cpt_dir,cpt_avg,cpt_set)=squeeze(Dcm(:,:,cpt_slc,cpt_b,cpt_dir,cpt_avg)).*exp(-T2Map(:,:,cpt_slc)/enum.TE(1));
                end
             end
           end
           waitbar(cpt_slc/size(enum.slc,2),h);
         end
    end
    close(h);    

end