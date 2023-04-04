function [fMap,VMap, DMap, S0Map, ErrorMap]= AIVIM_REDUCE_KM(Dcm, DcmB0, Bval,Aval,Db)

 warning off;
 
 fMap=[];
 VMap=[];
 DMap=[];
 S0Map=[];
 ErrorMap=[];
 
 disp('fIVIM fit') 
   

    
    
  for cpt=1:1:8
      
    if(cpt==1) % all
        b1=1;
        b2=7;
        a1=1;
        a2=5;
    elseif(cpt==2) % low a
        b1=1;
        b2=7;
        a1=1;
        a2=2;
    elseif(cpt==3) % high a
        b1=1;
        b2=7;
        a1=4;
        a2=5;
    elseif(cpt==4)% low b
        b1=1;
        b2=5;
        a1=1;
        a2=5;
    elseif(cpt==5) % mid b
        b1=2;
        b2=6;
        a1=1;
        a2=5;
    elseif(cpt==6) % high b
        b1=3;
        b2=7;
        a1=1;
        a2=5;   
    elseif(cpt==7) % IVIM 
        b1=1;
        b2=7;
        a1=5;
        a2=5;       
    else           % AIVIM
        b1=1;
        b2=7;
        a1=1;
        a2=5;
    end
      
    h = waitbar(0,['fIVIM fit num' num2str(cpt)]);  
    for y=2:1:(size(Dcm,1)-1)                          
        for x=2:1:(size(Dcm,2)-1)
            if DcmB0(y,x)>0
                
                Dat= squeeze( mean(mean(Dcm(y-1:y+1,x-1:x+1,1:a1:end,b1:b2))) );
                
                A1= reshape(Aval(a1:a2,b1:b2)',[(a2-a1+1)*(b2-b1+1) 1]);
                
                B1= repmat(Bval(b1:b2), 1, (a2-a1+1));
                
                D1= reshape(Dat' ,[(a2-a1+1)*(b2-b1+1) 1]);
                
                [f, V, D, ResM0]=  fitAFlowIVIM(D1', A1',B1,Db);
                
                S0Map(y,x,cpt)=0;
                fMap(y,x,cpt)=f;
                VMap(y,x,cpt)=V;   
                DMap(y,x,cpt)=D;
                            
                ErrorMap(y,x,cpt)=mean(ResM0);     
            else
                S0Map(y,x,cpt)=0;
                fMap(y,x,cpt)=0;
                VMap(y,x,cpt)=0;   
                DMap(y,x,cpt)=0;
                
                ErrorMap(y,x,cpt)=zeros(1,1);  
            end
        end
        waitbar(y/size(Dcm,1),h);
    end
    close(h); 
  end
  
  
  
  warning on;
    
end
       