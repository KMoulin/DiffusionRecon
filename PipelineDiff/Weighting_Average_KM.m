function [Dcm2 DcmB02 enum2]= Weighting_Average_KM(Dcm, DcmB0, enum, ErrorMap)

    enum2=enum; 
    Dcm2=[];
    DcmB02=[];
    disp('Weighting Average data') 
    h = waitbar(0,'Weighting Average data...');
    for cpt_slc=1:1:size(enum.slc,2)
       for cpt_b=1:1:size(enum.b,2)
           for cpt_dir=1:1:enum.dataset(cpt_b).dirNum
               tmpDcmData=[];
               if(cpt_b==1)
                 tmpDcmData=zeros(size(DcmB0(:,:,cpt_slc,1,:,:)));
               else
                 tmpDcmData=zeros(size(Dcm(:,:,cpt_slc,1,1,1)));  
               end
               
                 if(cpt_b==1)
                     tmpDcmData=sum(squeeze(DcmB0(:,:,cpt_slc,cpt_b,cpt_dir,:)),3);
                 else
                     tmpDcmData=sum(squeeze(Dcm(:,:,cpt_slc,(cpt_b-1),cpt_dir,:)).*squeeze(ErrorMap(:,:,cpt_slc,(cpt_b-1),cpt_dir,:)),3);
                 end
                 if(cpt_b==1)
                     DcmB02(:,:,cpt_slc,cpt_b,cpt_dir,1)=tmpDcmData/enum.dataset(cpt_b).dir(cpt_dir).avgNum;
                 else
                     Dcm2(:,:,cpt_slc,(cpt_b-1),cpt_dir,1)=tmpDcmData./squeeze(sum(ErrorMap(:,:,cpt_slc,(cpt_b-1),cpt_dir,:),6));
                 end
               enum2.dataset(cpt_b).dir(cpt_dir).avgNum=1;
           end
       end
       waitbar(cpt_slc/size(enum.slc,2),h);
    end
    close(h);    

end