function [fMap,VMap, DMap, S0Map, ErrorMap, IndexMap]= Rand_AIVIM_KM(Dcm, fSpace,Db, NumRand,NumTab,NumPoint)
    warning off
    fMap=[];
    VMap=[];
    DMap=[];
    S0Map=[];
    ErrorMap=[];
    
    lb =   [0  0  0];    
    ub =   [1   0.1   100];     
    val0 = [0.15   0.001     4];   
    options=optimset('lsqcurvefit');
    options.Algorithm='levenberg-marquardt';
    options.TolFun=1e-8;
    options.TolX=1e-8;
    options.Display='off';
  
    IVIM_M1_relax = @(x,xdata)(  (1- x(1)).*exp( -x(2).*xdata(1,:) )   +     x(1).*exp( -Db.*xdata(1,:) ).* exp( -xdata(2,:).*x(3).*x(3)/2 ) );
    
    Alpha=fSpace(:,1)';
    Bval=fSpace(:,2)';
    
    disp('IVIM fit') 
    h = waitbar(0,['IVIM fit ' num2str(1)]);
    for t=1:1:NumRand
        
       
        index=NumTab(t,1:NumPoint);
        t
        for y=2:1:(size(Dcm,1)-1)                          
            for x=2:1:(size(Dcm,2)-1)
                    Dat = squeeze( mean(mean(Dcm(y-1:y+1,x-1:x+1,:,:))) );
                    %Dat= squeeze( ((Dcm(y,x,:,:))) );
                    Dat=reshape(squeeze(Dat)',size(Dcm,3)*size(Dcm,4),1);
                    
                    [val,resnorm,residual,exitflag,output] = lsqcurvefit(IVIM_M1_relax, val0, [Bval(index); Alpha(index)], Dat(index)', lb, ub, options);
                    
                    fMap(y,x,t)=val(1);
                    VMap(y,x,t)=val(3);
                    DMap(y,x,t)=val(2);
                    S0Map(y,x,t)=0;
                    ErrorMap(y,x,t,:)=residual./Dat(index)';
                    
                  
                    IndexMap(t,1,:)=ones(length(Alpha),1)*9999;
                    IndexMap(t,1,index)=Alpha(index);
                    IndexMap(t,2,:)=ones(length(Bval),1)*9999;
                    IndexMap(t,2,index)=Bval (index);
                    
                    
            end
            waitbar(y/size(Dcm,1),h);
        end
        
    end
    close(h);
    warning on;
end
       