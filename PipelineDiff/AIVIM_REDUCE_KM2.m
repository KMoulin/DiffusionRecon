function [fMap,VMap, DMap, S0Map, ErrorMap]= AIVIM_REDUCE_KM2(Dcm, DcmB0, Bval,Aval,Db)

 warning off;
 
 fMap=[];
 VMap=[];
 DMap=[];
 S0Map=[];
 ErrorMap=[];
 
 disp('fIVIM fit') 
   
    h = waitbar(0,['fIVIM fit num' num2str(0)]);  
    for y=2:1:(size(Dcm,1)-1)                          
        for x=2:1:(size(Dcm,2)-1)
            if DcmB0(y,x)>0
                
                Dat= squeeze( mean(mean(Dcm(y-1:y+1,x-1:x+1,end,:))) );
                
                A1= reshape(Aval(end,:)',[7 1]);
                
                B1= repmat(Bval, 1, 1);
                
                D1= reshape(Dat' ,[7 1]);
                
                [f, V, D, ResM0]=  fitAFlowIVIM(D1', A1',B1,Db);
                
                S0Map(y,x,1)=0;
                fMap(y,x,1)=f;
                VMap(y,x,1)=V;   
                DMap(y,x,1)=D;
                            
                ErrorMap(y,x,1)=mean(ResM0);     
            else
                S0Map(y,x,1)=0;
                fMap(y,x,1)=0;
                VMap(y,x,1)=0;   
                DMap(y,x,1)=0;
                
                ErrorMap(y,x,1)=zeros(1,1);  
            end
        end
        waitbar(y/size(Dcm,1),h);
    end
    close(h); 
    
    
    h = waitbar(0,['fIVIM fit num' num2str(1)]);  
    for y=2:1:(size(Dcm,1)-1)                          
        for x=2:1:(size(Dcm,2)-1)
            if DcmB0(y,x)>0
                
                Dat= squeeze( mean(mean(Dcm(y-1:y+1,x-1:x+1,1:4:end,:))) );
                
                A1= reshape(Aval(1:4:end,:)',[14 1]);
                
                B1= repmat(Bval, 1, 2);
                
                D1= reshape(Dat' ,[14 1]);
                
                [f, V, D, ResM0]=  fitAFlowIVIM(D1', A1',B1,Db);
                
                S0Map(y,x,2)=0;
                fMap(y,x,2)=f;
                VMap(y,x,2)=V;   
                DMap(y,x,2)=D;                          
                ErrorMap(y,x,2)=mean(ResM0);     
                
            else
                S0Map(y,x,2)=0;
                fMap(y,x,2)=0;
                VMap(y,x,2)=0;   
                DMap(y,x,2)=0;
                ErrorMap(y,x,2)=zeros(1,1);  
                
            end
        end
        waitbar(y/size(Dcm,1),h);
    end
    close(h); 
  warning on;
    
end
       