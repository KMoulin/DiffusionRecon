function [Eval]= Eval_KM2(Dcm, DcmB0, enum)
 
    disp('Eval Map') 
    Eval=[];
    h = waitbar(0,'Error Map...');
    for cpt_slc=1:1:size(enum.slc,2)
       for cpt_b=1:1:size(enum.b,2)
           k=1;
           
           for cpt_dir=1:1:(enum.dataset(cpt_b).dirNum)
               if(cpt_b~=1)
                     Ref=max(Dcm(:,:,cpt_slc,(cpt_b-1),cpt_dir,:),[],6);
                     Ref(DcmB0(:,:,cpt_slc)==0)=0;
                     RefD=repmat(Ref,1,1,5);
                     xEval=linspace(0,Extract_KM(Ref),100);
                     for cpt_eval=1:1:100
                        Tmp=squeeze(Dcm(:,:,cpt_slc,(cpt_b-1),cpt_dir,:));
                        Tmp(Tmp<=RefD.*cpt_eval/100)=0;
                        TmpAvg=sum(Tmp,3)./sum(Tmp ~= 0,3); 
                        TmpADC =ADCMap_KM(TmpAvg,squeeze(DcmB0(:,:,cpt_slc)),-enum.b(cpt_b));
                        Eval(cpt_slc,(cpt_b-1),cpt_dir,cpt_eval)=Extract_KM(TmpADC) ;                  
                     end
               end
           end        
       end
       waitbar(cpt_slc/size(enum.slc,2),h); 
    end
    close(h);
    

end