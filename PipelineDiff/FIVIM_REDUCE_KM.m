function [fMap,VMap,S0Map, ErrorMap]= FIVIM_REDUCE_KM(Dcm, DcmB0, Bval,Aval,Db)

    fMap=[];
    VMap=[];
    S0Map=[];
    ErrorMap=[];
    disp('fIVIM fit') 
   
    h = waitbar(0,['fIVIM fit slc']);
    
    
    for b=(size(Dcm,3)):-1:3
    for y=2:1:(size(Dcm,1)-1)                          
        for x=2:1:(size(Dcm,2)-1)
         
            if DcmB0(y,x)>0
                Dat=squeeze(mean(mean(Dcm(y-1:y+1,x-1:x+1,1:b,:))));
                a1=reshape(Aval(:,2:(b+1))',[4*(b) 1]);
                b1=[Bval(2:b+1) Bval(2:b+1) Bval(2:b+1) Bval(2:b+1)];
                d1=reshape(Dat ,[4*(b) 1]);
                [f, V, ResM0]=  fitFFlowIVIM(d1', a1',b1,Db);
                S0Map(y,x,b)=0;
                fMap(y,x,b)=f;
                VMap(y,x,b)=V;     
               % ErrorMap(y,x,b,:)=ResM0;     
            else
                S0Map(y,x,b)=0;
                fMap(y,x,b)=0;
                VMap(y,x,b)=0;
               % ErrorMap(y,x,b,:)=zeros(1,1,size(Dcm,4));  
       
           end
        end
        waitbar(y/size(Dcm,1),h);
    end
    end
    close(h); 
end
       