function [Dcm2 DcmB02]= Filter_KM(Dcm, DcmB0, enum,h)

    enum2=enum; 
    Dcm2=[];
    DcmB02=[];
    disp('tMIP data') 
    h = waitbar(0,'tMIP data...');
    for cpt_slc=1:1:size(enum.slc,2)
       for cpt_b=1:1:size(enum.b,2)
           for cpt_dir=1:1:enum.dataset(cpt_b).dirNum;           
             if(cpt_b==1)DcmB02(:,:,cpt_slc,cpt_b,cpt_dir,1)=max(DcmB0(:,:,cpt_slc,cpt_b,cpt_dir,:),[],6);
             else Dcm2(:,:,cpt_slc,(cpt_b-1),cpt_dir,1)=max(Dcm(:,:,cpt_slc,(cpt_b-1),cpt_dir,:),[],6);
             end
             enum2.dataset(cpt_b).dir(cpt_dir).avgNum=1;
           end
       end
       waitbar(cpt_slc/size(enum.slc,2),h);
    end
    close(h);    

end