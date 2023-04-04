function [Dcm2 DcmB02 enum2]= Trace_and_Reject_KM(Dcm, DcmB0, enum)

% Perform averaging on the DWI and nDWI matrices
% 
% Update the number of averages in enum
%
%
% SYNTAX:  [Dcm2 DcmB02 enum2]= Average_and_Reject_KM(Dcm, DcmB0, enum);
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
% Kevin Moulin 08.31.2017
% Kevin.Moulin.26@gmail.com
% Ennis Lab @ UCLA; http://mrrl.ucla.edu

    enum2=enum; 
    Dcm2=[];
    DcmB02=[];
    disp('Average and reject data') 
    h = waitbar(0,'Average and reject data...');
    for cpt_slc=1:1:size(enum.slc,2)
       for cpt_b=1:1:size(enum.b,2)
           if(cpt_b==1)
                for cpt_dir=1:1:enum.dataset(cpt_b).dirNum;
                    DcmB02(:,:,cpt_slc,cpt_b,cpt_dir,1)=mean(DcmB0(:,:,cpt_slc,cpt_b,cpt_dir,:),6);
                end
                enum2.dataset(cpt_b).dir(cpt_dir).avgNum=1;
           else          
                        tmp_Dcm=squeeze(Dcm(:,:,cpt_slc,(cpt_b-1),:,:));
                        tmp_DcmB0=squeeze(DcmB02(:,:,cpt_slc,1,1,:));
                        tmpADC=log(tmp_Dcm./repmat(tmp_DcmB0,1,1,enum2.dataset(cpt_b).dirNum))/-enum.b(cpt_b);
                        tmpADC(tmpADC>0.003)=0;    
                        tmpAvg=tmpADC;
                        tmpAvg(tmpAvg>0.000001)=1;
                        Mask=mean(tmpAvg,3);
                        Dcm2(:,:,cpt_slc,(cpt_b-1),1,1)=mean(squeeze(Dcm(:,:,cpt_slc,(cpt_b-1),:,:)).*tmpAvg,3)./Mask;
                                 
           end
       end
       waitbar(cpt_slc/size(enum.slc,2),h);
    end
    for cpt_slc=1:1:size(enum.slc,2)
       for cpt_b=1:1:size(enum.b,2)
           enum2.dataset(cpt_b).dirNum=1;
       end
    end
    close(h);    

end