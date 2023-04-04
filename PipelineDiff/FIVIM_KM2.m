function [fMap,VMap,S0Map, ErrorMap]= FIVIM_KM2(Dcm, Mask, fSpace,Db)

    fMap=[];
    VMap=[];
    S0Map=[];
    ErrorMap=[];
    disp('fIVIM fit') 
   
    h = waitbar(0,['fIVIM fit slc']);
    
    
    warning off;
    for y=1:1:(size(Dcm,1))                          
        for x=1:1:(size(Dcm,2))
            if Mask(y,x)>0
                Dat=squeeze(((Dcm(y,x,:))))';
                [f, V, ResM0]=  fitFFlowIVIM(Dat,fSpace(:,1)',fSpace(:,2)',Db); % Alpha, Bval
                S0Map(y,x)=0;
                fMap(y,x)=f;
                VMap(y,x)=V;     
                ErrorMap(y,x,:)=0;     
            else
                S0Map(y,x)=0;
                fMap(y,x)=0;
                VMap(y,x)=0;
                ErrorMap(y,x,:)=zeros(1,1,size(Dcm,3));  
            end
        end
        waitbar(y/size(Dcm,1),h);
    end
    close(h); 
    warning on;
end
       