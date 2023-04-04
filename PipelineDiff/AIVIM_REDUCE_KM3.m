function [fMap,VMap, DMap, S0Map, ErrorMap]= AIVIM_REDUCE_KM3(Dcm, DcmB0, Bval,Aval,Db)

 warning off;
 
 fMap=[];
 VMap=[];
 DMap=[];
 S0Map=[];
 ErrorMap=[];
 
 disp('fIVIM fit') 
   

    
      A1= reshape(Aval',[35 1]);             
      B1= repmat(Bval, 1,5);
  
      ErrorMap  = zeros(size(Dcm,1) , size(Dcm,2) ,5, 7);
      for cpt=35:-34:4 % -1

        h = waitbar(0,['fIVIM fit num' num2str(cpt)]);  
        for y=2:1:(size(Dcm,1)-1)                          
            for x=2:1:(size(Dcm,2)-1)
                if DcmB0(y,x)>0

                    Dat= squeeze( mean(mean(Dcm(y-1:y+1,x-1:x+1,:,:))) );
                    D1= reshape(Dat' ,[35, 1]);
                    D1(B1==0)=1;
                    
                    
                    [f, V, D, ResM0]=  fitAFlowIVIM(D1', A1',B1,Db);

                    S0Map(y,x,cpt)=0;
                    fMap(y,x,cpt)=f;
                    VMap(y,x,cpt)=V;   
                    DMap(y,x,cpt)=D;

                    ErrorMap(y,x,:,:)=reshape((ResM0),[5, 7]);    
                    if(cpt==35)
                        RefErrorMap(y,x,:)=(ResM0);
                    end
                else
                    S0Map(y,x,cpt)=0;
                    fMap(y,x,cpt)=0;
                    VMap(y,x,cpt)=0;   
                    DMap(y,x,cpt)=0; 
                end
            end
            waitbar(y/size(Dcm,1),h);
        end
        close(h);
%         if(cpt>4)
%             %[M I]=max(abs(squeeze(mean(mean(ErrorMap(:,:,cpt,:))))));
%             
%             [M I]=max(squeeze(mean(mean(abs(RefErrorMap)))));
%             RefErrorMap(:,:,I)=zeros(size(RefErrorMap,1),size(RefErrorMap,2));
%             B1(I)=0;
%             A1(I)=0;
%         end
%         
        
      end



      warning on;

end
       