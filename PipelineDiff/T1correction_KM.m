function Dcm2= T1correction_KM(Dcm,T1Map, Library,Index,enum)


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

 
  
    disp('T1 Correction') 
    h = waitbar(0,'T1 Correction...');
    for cpt_set=1:1:1%enum.nset
         for cpt_slc=1:1:enum.datasize(cpt_set).slc
            for cpt_b=1:1:enum.datasize(cpt_set).b     
              for cpt_dir=1:1: enum.dataset(cpt_set).slc(cpt_slc).b(cpt_b).nb_dir           
                 for cpt_avg=1:1: enum.dataset(cpt_set).slc(cpt_slc).b(cpt_b).dir(cpt_dir).nb_avg   
                     for y=1:1:size(Dcm,1)
                         for x=1:1:size(Dcm,2)
                             if T1Map(y,x,cpt_slc)~=0
                                Dcm2(y,x,cpt_slc,cpt_b,cpt_dir,cpt_avg,cpt_set)=Dcm(y,x,cpt_slc,cpt_b,cpt_dir,cpt_avg)./squeeze(abs(Library(cpt_slc,cpt_b,cpt_dir,cpt_avg,find(Index(:,1)==T1Map(y,x,cpt_slc)))));
                             else
                                 Dcm2(y,x,cpt_slc,cpt_b,cpt_dir,cpt_avg,cpt_set)=0;
                             end
                         end
                     end
                end
             end
           end
           waitbar(cpt_slc/size(enum.slc,2),h);
         end
    end
    close(h);    


%     Dcm=squeeze(Dcm);  
%     for y=1:1:size(Dcm,1)
%         for x=1:1:size(Dcm,2)
%             for z=1:1:size(Dcm,3)
%                     if T1Map(y,x,z)~=0
%                         Dcm2(y,x,z,:)=squeeze(Dcm(y,x,z,:))./squeeze(abs(Library(z,1,1,:,find(Index(:,1)==T1Map(y,x,z)))));
%                     else
%                         Dcm2(y,x,z,:)=zeros(size(Dcm,4),1);
%                     end
%             end
%         end
%     end

end