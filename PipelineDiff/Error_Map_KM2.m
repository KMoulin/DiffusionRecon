function [Dcm2 ErrorMap2]= Error_Map_KM2(Dcm, DcmB0, enum)
 
    Dcm2=[];
    disp('Error Map') 
    ErrorMap=zeros(size(Dcm,1),size(Dcm,2),size(enum.slc,2),enum.dataset(2).dirNum,enum.dataset(2).dir(1).avgNum,3);
    ErrorMap2=[];
    for cpt_slc=1:1:size(enum.slc,2)
       for cpt_b=1:1:size(enum.b,2)
           k=1;
           h = waitbar(0,'Error Map...');
           for cpt_dir=1:1:(enum.dataset(cpt_b).dirNum)
               if(cpt_b~=1)
                    Ref=max(Dcm(:,:,cpt_slc,(cpt_b-1),cpt_dir,:),[],6);
               end
             for cpt_avg=1:1:(enum.dataset(cpt_b).dir(cpt_dir).avgNum)
                if(cpt_b~=1)

                         Tmp=squeeze(Dcm(:,:,cpt_slc,(cpt_b-1),cpt_dir,cpt_avg));
                         ErrorMap2(:,:,cpt_slc,cpt_dir,cpt_avg)=Tmp./Ref ;
                         k=k+1;
                end
             end 
            waitbar(cpt_dir/(enum.dataset(cpt_b).dirNum),h); 
            end
           close(h);
       end
       
    end
    
   
%      for cpt_slc=1:1:size(enum.slc,2)
%          
%        DxMin=0.0005;%min(Dcm2(:,:,cpt_slc,:,1),[],4);
%        DyMin=0.0005;%min(Dcm2(:,:,cpt_slc,:,2),[],4);
%        DzMin=0.0005;%min(Dcm2(:,:,cpt_slc,:,3),[],4);
% 
%        for cpt_b=1:1:size(enum.b,2)      
%            for cpt_dir=1:1:(enum.dataset(cpt_b).dirNum);                              
%                  for cpt_avg=1:1:(enum.dataset(cpt_b).dir(cpt_dir).avgNum)
%                    if(cpt_b~=1)    
%                     
%                     ErrorMap(:,:,cpt_slc,cpt_dir,cpt_avg,1)=(1-DxMin./ErrorMap(:,:,cpt_slc,cpt_dir,cpt_avg,1));
%                     ErrorMap(:,:,cpt_slc,cpt_dir,cpt_avg,2)=(1-DyMin./ErrorMap(:,:,cpt_slc,cpt_dir,cpt_avg,2));
%                     ErrorMap(:,:,cpt_slc,cpt_dir,cpt_avg,3)=(1-DzMin./ErrorMap(:,:,cpt_slc,cpt_dir,cpt_avg,3));   
%                        
%                     ErrorMap2(:,:,cpt_slc,cpt_dir,cpt_avg)=(ErrorMap(:,:,cpt_slc,cpt_dir,cpt_avg,1)+ErrorMap(:,:,cpt_slc,cpt_dir,cpt_avg,2)+ErrorMap(:,:,cpt_slc,cpt_dir,cpt_avg,3))/3 ;
%                    end  
%                  end
%            end
%        end
%      end
%      
    

end