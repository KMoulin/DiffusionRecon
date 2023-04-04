function [fMap,VMap, DMap, S0Map, ErrorMap]= AIVIM_KM(Dcm, Mask, fSpace,Db)

    fMap=[];
    VMap=[];
    DMap=[];
    S0Map=[];
    ErrorMap=[];
    disp('AIVIM fit') 
   
    h = waitbar(0,['AIVIM fit']);
    
    for y=2:1:(size(Dcm,1)-1)                          
        for x=2:1:(size(Dcm,2)-1)
            if Mask(y,x)>0
               % Dat= squeeze( mean(mean(Dcm(y-1:y+1,x-1:x+1,:,:))) );
                Dat= squeeze( ((Dcm(y,x,:,:))) );
                Dat=reshape(squeeze(Dat)',size(Dcm,3)*size(Dcm,4),1);
                [f, V, D, ResM0]=  fitAFlowIVIM(Dat',(fSpace(:,1)').*(fSpace(:,1)'),fSpace(:,2)',Db); % Alpha, Bval
                fMap(y,x)=f;
                VMap(y,x)=V;
                DMap(y,x)=D;
                S0Map(y,x)=0;
                ErrorMap(y,x,:)=ResM0;     
            else
                fMap(y,x)=0;
                VMap(y,x)=0;
                DMap(y,x)=0;
                S0Map(y,x)=0;
                ErrorMap(y,x,:)=zeros(1,1,size(Dcm,3));  
            end
        end
        waitbar(y/size(Dcm,1),h);
    end
    close(h); 
end
       