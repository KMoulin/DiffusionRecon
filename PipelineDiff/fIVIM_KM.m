function [fMap,VMap,DMap,S0Map, ErrorMap]= fIVIM_KM(Dcm, Mask, b_vect,a_vect)

    fMap=[];
    VMap=[];
    DMap=[];
    S0Map=[];
    ErrorMap=[];
    disp('fIVIM fit') 
   
    h = waitbar(0,'fIVIM fit slc');
    
    
    warning off;
    z=1;
    d=1;
   % for d=1:1:(size(Dcm,5))
        for y=2:1:(size(Dcm,1)-1)                          
            for x=2:1:(size(Dcm,2)-1)
                %for z=1:1:(size(Dcm,3))

                        if Mask(y,x,z)>0
                            %Dat=squeeze(((Dcm(y,x,z,:,d))))';
                            Dat= squeeze( mean(mean(Dcm(y-1:y+1,x-1:x+1,:,:))) );
                           [f, V, D, ResM0]=  fitAFlowIVIM(Dat,a_vect(d,2:end),b_vect,1.75e-3);
                          %  [f, D, Dstar, ResM0]=  fit2pointLB(Dat,b_vect) ;
                            S0Map(y,x,z,d)=0;
                            fMap(y,x,z,d)=f;
                            VMap(y,x,z,d)=V;
                            DMap(y,x,z,d)=D;  
                            ErrorMap(y,x,z,d,:)=0;     
                        else
                            S0Map(y,x,z,d)=0;
                            fMap(y,x,z,d)=0;
                            VMap(y,x,z,d)=0;
                            DMap(y,x,z,d)=0;
                            ErrorMap(y,x,z,d,:)=zeros(1,1,size(Dcm,4));  
                        end        
               % end
            end
            waitbar(y/size(Dcm,1),h);
        end
   % end
    close(h); 
    warning on;
end
       