function [Dcm2]= RigidRegistration_KM2(Dcm, enum)

% Register the DWI and nDWI matrices per SLICE based on a rigid method.
% Registration reference is automatically chosen to be the image with the
% maximum signal intensity (usually the nDWI b=0 s/mm² images)
%
% SYNTAX:  [Dcm2]= RigidRegistration_KM(Dcm, enum);
%  
%
% INPUTS:   Dcm - DWI image matrix
%                 [y x slices b-values directions averages]        
%           
%           enum - Structure which contains information about the dataset 
%          
% OUTPUTS:  Dcm2 - DWI image matrix 
%                 [y x slices b-values directions averages]
%                      
%
% Kevin Moulin 08.14.2017
% Kevin.Moulin.26@gmail.com
% Ennis Lab @ UCLA; http://mrrl.ucla.edu
   
    Dcm2=[];
    [optimizer, metric] = imregconfig('multimodal');
    disp('Rigid registration before') 
    h = waitbar(0,'Rigid registration before...');
    for cpt_set=1:1:enum.nset
        for cpt_slc=1:1:enum.datasize(cpt_set).slc              
            for cpt_b=1:1:enum.datasize(cpt_set).b     
               for cpt_dir=1:1: enum.dataset(cpt_set).slc(cpt_slc).b(cpt_b).nb_dir  
                  ref=zeros(size(Dcm,1),size(Dcm,2)); 
                  for cpt_avg=1:1: enum.dataset(cpt_set).slc(cpt_slc).b(cpt_b).dir(cpt_dir).nb_avg  
                        if (mean(mean(ref))<mean(mean(Dcm(:,:,cpt_slc,cpt_b,cpt_dir,cpt_avg,cpt_set))))
                            ref=Dcm(:,:,cpt_slc,cpt_b,cpt_dir,cpt_avg,cpt_set);
                        end                          
                  end
                  for cpt_avg=1:1: enum.dataset(cpt_set).slc(cpt_slc).b(cpt_b).dir(cpt_dir).nb_avg  
                            Dcm2(:,:,cpt_slc,cpt_b,cpt_dir,cpt_avg,cpt_set)=  imregister(Dcm(:,:,cpt_slc,cpt_b,cpt_dir,cpt_avg,cpt_set), ref, 'affine' , optimizer, metric); 
                  end
               end
            end
%             for cpt_b=1:1:enum.datasize(cpt_set).b     
%                for cpt_dir=1:1: enum.dataset(cpt_set).slc(cpt_slc).b(cpt_b).nb_dir           
%                   for cpt_avg=1:1: enum.dataset(cpt_set).slc(cpt_slc).b(cpt_b).dir(cpt_dir).nb_avg  
%                             Dcm2(:,:,cpt_slc,cpt_b,cpt_dir,cpt_avg,cpt_set)=  imregister(Dcm(:,:,cpt_slc,cpt_b,cpt_dir,cpt_avg,cpt_set), ref, 'affine', optimizer, metric); 
%                   end           
%                end
%             end
           waitbar(cpt_slc/size(enum.slc,2),h);
        end
    end
   
    close(h);  

end