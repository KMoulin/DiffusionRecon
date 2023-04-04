function [T1map T2map Offmap]= T1T2fit_KM(Dcm, enum, Library, Index)

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

    Dcm=squeeze(Dcm);
    %Dcm2=squeeze(Dcm./repmat(max(Dcm(:,:,:,:),[],4),1,1,1,size(Dcm,4)));
     Dcm2=squeeze(Dcm./repmat(Dcm(:,:,:,6),1,1,1,size(Dcm,4)));
    Library2=Library;
   % Dcm2=Dcm2./repmat(mean(Dcm2(:,:,:,end-enum.nRestore:end),4),1,1,1,14);
    
   % Library2=Library(:,:,:,(end+1-enum.nRestore:end),:);%end+1-enum.nRestore
  
   Library2=Library2./repmat(Library2(:,:,:,6,:),1,1,1,size(Library2,4),1);
    
  %  Library2=Library./repmat(max(Library(:,:,:,:,:),[],4),1,1,1,size(Library,4),1);
    T1map=[];
    T2map=[];
    Offmap=[];
    disp('T1 fit') 
    h = waitbar(0,'T1 fit...');
    for cpt_slc=1:1:size(Dcm,3)
        tic
        for cpt_y=1:1:size(Dcm,1)
            for cpt_x=1:1:size(Dcm,2)        
                if (squeeze(Dcm2(cpt_y,cpt_x,cpt_slc,1))>0)
                    VectResult=sqrt(mean(   abs(squeeze((Library2(cpt_slc,1,1,6:end,:)))) - squeeze(Dcm2(cpt_y,cpt_x,cpt_slc,6:end))   ).^2);
                    [ValMin IdxMin]=min(VectResult);
                    T1map(cpt_y,cpt_x,cpt_slc)=Index(IdxMin,1);
                    T2map(cpt_y,cpt_x,cpt_slc)=Index(IdxMin,2);
                    Offmap(cpt_y,cpt_x,cpt_slc)=Index(IdxMin,3);
                else
                    T1map(cpt_y,cpt_x,cpt_slc)=0;
                    T2map(cpt_y,cpt_x,cpt_slc)=0;
                    Offmap(cpt_y,cpt_x,cpt_slc)=0;
                end
            end
        end
        waitbar(cpt_slc/size(Dcm,3),h);
        toc
    end
    close(h);




%% 
cpt_x=55;
cpt_y=74;

cpt_x=55;
cpt_y=63;
for cpt_slc=1:1:5
VectResult=sqrt(mean(   squeeze(abs(Library2(cpt_slc,1,1,6:end,:))) - squeeze(Dcm2(cpt_y,cpt_x,cpt_slc,6:end))   ).^2);
[ValMin IdxMin]=min(VectResult);
Index(IdxMin,1);
Index(IdxMin,2);
figure,plot(squeeze(abs(Library2(cpt_slc,1,1,:,IdxMin)))'),hold on,plot(squeeze(Dcm2(cpt_y,cpt_x,cpt_slc,:))')
end

end