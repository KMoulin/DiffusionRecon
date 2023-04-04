function [fMap,VMap,S0Map, ErrorMap]= FIVIM_ALPHA_KM(Dcm, DcmB0, Bval,Aval,Db)

    fMap=[];
    VMap=[];
    S0Map=[];
    ErrorMap=[];
    disp('fIVIM fit') 
   
    h = waitbar(0,['fIVIM fit slc']);
    
    
    
    for y=2:1:(size(Dcm,1)-1)                          
        for x=2:1:(size(Dcm,2)-1)
           for b=1:1:size(Dcm,3) 
            if DcmB0(y,x)>0
                Dat=squeeze(mean(mean(Dcm(y-1:y+1,x-1:x+1,b,:))));
                [f, V, ResM0]=  fitFFlowAlphaIVIM(Dat, Aval,Bval(b+1),Db);
                S0Map(y,x,b)=0;
                fMap(y,x,b)=f;
                VMap(y,x,b)=V;     
                ErrorMap(y,x,b,:)=ResM0;     
            else
                S0Map(y,x,b)=0;
                fMap(y,x,b)=0;
                VMap(y,x,b)=0;
                ErrorMap(y,x,b,:)=zeros(1,1,size(Dcm,4));  
            end
           end
        end
        waitbar(y/size(Dcm,1),h);
    end
    close(h); 
end
       