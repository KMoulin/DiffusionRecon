function [fMap,VMap, DMap, S0Map,ErrorMap]= AIVIMB0_KM(Dcm,DcmB0, Bval,Aval,Db)

    fMap=[];
    VMap=[];
    DMap=[];
    S0Map=[];
    ErrorMap=[];
    disp('fIVIM fit') 
   
    h = waitbar(0,['fIVIM fit slc']);
    
    
    
    for y=1:1:(size(Dcm,1))                          
        for x=1:1:(size(Dcm,2))
           if DcmB0(y,x,1)>0
            [D, S0, ResM0]=  fitADCB0(squeeze(Dcm(y,x,:,:))',Bval) ;
            %[f, V, D, S0, ResM0]=  fitAFlowIVIMB0(squeeze(Dcm(y,x,:,:))',Aval,Bval,Db);
           % fMap(y,x)=f;
           % VMap(y,x)=V;
            DMap(y,x)=D;
            S0Map(y,x)=S0;
            ErrorMap(y,x,:)=ResM0;   
           
           else
            DMap(y,x)=0;
            S0Map(y,x)=0;
            ErrorMap(y,x,1:size(Dcm,4))=0;      
           end
        end
        waitbar(y/size(Dcm,1),h);
    end
    close(h); 
end
       