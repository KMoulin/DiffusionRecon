function [Dcm2 DcmB02]= RigidRegistration_IVIM_KM(Dcm, DcmB0, enum)

    Dcm2=[];
    DcmB02=[];
    DcmM02=[];
    DcmB0M02=[];
    [optimizer, metric] = imregconfig('multimodal');
    disp('Rigid registration') 
    h = waitbar(0,'Rigid registration...');
    for cpt_slc=1:1:size(enum.slc,2)
       ref(:,:,cpt_slc)=zeros(size(DcmB0,1),size(DcmB0,2));
       for cpt_b=1:1:size(enum.b,2)
           for cpt_dir=1:1:enum.dataset(cpt_b).dirNum;
               for cpt_avg=1:1:enum.dataset(cpt_b).dir(cpt_dir).avgNum
                 if(cpt_b==1)                   
                         if (mean(mean(ref(:,:,cpt_slc)))<mean(mean(DcmB0(:,:,cpt_slc,cpt_b,cpt_dir,cpt_avg))))
                             ref(:,:,cpt_slc)=DcmB0(:,:,cpt_slc,cpt_b,cpt_dir,cpt_avg);
                         end
                 else
                        if (mean(mean(ref(:,:,cpt_slc)))<mean(mean(Dcm(:,:,cpt_slc,cpt_b-1,cpt_dir,cpt_avg))))
                            ref(:,:,cpt_slc)=Dcm(:,:,cpt_slc,cpt_b-1,cpt_dir,cpt_avg);
                        end
                 end           
               end           
           end
       end
      
    end

    for cpt_slc=1:1:size(enum.slc,2)
       for cpt_b=1:1:size(enum.b,2)
           for cpt_dir=1:1:enum.dataset(cpt_b).dirNum;
               for cpt_avg=1:1:enum.dataset(cpt_b).dir(cpt_dir).avgNum
                 if(cpt_b==1)                   
                        DcmB02(:,:,cpt_slc,cpt_b,cpt_dir,cpt_avg)=  imregister(DcmB0(:,:,cpt_slc,cpt_b,cpt_dir,cpt_avg), ref(:,:,1), 'affine', optimizer, metric);
                 else
                        Dcm2(:,:,cpt_slc,cpt_b-1,cpt_dir,cpt_avg)=  imregister(Dcm(:,:,cpt_slc,cpt_b-1,cpt_dir,cpt_avg), ref(:,:,1), 'affine', optimizer, metric);
                 end           
               end           
           end
       end
       waitbar(cpt_slc/size(enum.slc,2),h);
    end
    close(h);  

end