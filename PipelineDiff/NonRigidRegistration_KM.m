function [Dcm2]= NonRigidRegistration_KM(Dcm, enum)

% Register the DWI and nDWI matrices per SLICE based on a nonrigid method.
% Registration reference is automatically chosen to be the image with the
% maximum signal intensity (usually the nDWI b=0 s/mm² images). Reject 
% images that are not properlly register and generate blank images instead
%
% SYNTAX:  [Dcm2]= NonRigidRegistration_KM(Dcm, enum);
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
    disp('NonRigid registration') 
    h = waitbar(0,'NonRigid registration...');
    for cpt_set=1:1:enum.nset
        for cpt_slc=1:1:enum.datasize(cpt_set).slc
            ref=zeros(size(Dcm,1),size(Dcm,2));   
            for cpt_b=1:1:enum.datasize(cpt_set).b     
               for cpt_dir=1:1: enum.dataset(cpt_set).slc(cpt_slc).b(cpt_b).nb_dir           
                  for cpt_avg=1:1: enum.dataset(cpt_set).slc(cpt_slc).b(cpt_b).dir(cpt_dir).nb_avg  
                        if (mean(mean(ref))<mean(mean(Dcm(:,:,cpt_slc,cpt_b,cpt_dir,cpt_avg,cpt_set))))
                            ref=Dcm(:,:,cpt_slc,cpt_b,cpt_dir,cpt_avg,cpt_set);
                        end                          
                  end           
               end
            end
            for cpt_b=1:1:enum.datasize(cpt_set).b     
               for cpt_dir=1:1: enum.dataset(cpt_set).slc(cpt_slc).b(cpt_b).nb_dir           
                  for cpt_avg=1:1: enum.dataset(cpt_set).slc(cpt_slc).b(cpt_b).dir(cpt_dir).nb_avg 
                        [Mp,eng] = DiffRegistrationLinear(ref,Dcm(:,:,cpt_slc,cpt_b,cpt_dir,cpt_avg,cpt_set));%%DiffRegistrationLinear
                        if eng<550
                            Dcm2(:,:,cpt_slc,cpt_b,cpt_dir,cpt_avg)=Mp;
                        else
                            Dcm2(:,:,cpt_slc,cpt_b,cpt_dir,cpt_avg)=zeros(size(Mp,1),size(Mp,2)); 
                           % cpt_avg_tmp=cpt_avg_tmp-11;
                        end     
                  end           
               end
            end
           waitbar(cpt_slc/size(enum.slc,2),h);
        end
    end
   
    close(h);  

end