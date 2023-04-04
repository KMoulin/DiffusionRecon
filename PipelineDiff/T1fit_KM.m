function [T1Map, T1starMap]= T1fit_KM(Dcm, enum ,TI)

    T1Map=[];
    T1starMap=[];
    AMap=[];
    BMap=[];
    resMap=[];
    
    disp('T1 fit') 
   
    
    warning off;
    for slc=1:1:(size(Dcm,3))
        tic
        h = waitbar(0,['T1 fit slc' num2str(slc) ]);
        for y=2:1:(size(Dcm,1)-1)                          
            for x=2:1:(size(Dcm,2)-1)
                        if Dcm(y,x,slc,1)>0
                            Dat=squeeze(((Dcm(y,x,slc,:))));
                           % Dat= squeeze( mean(mean(Dcm(y-1:y+1,x-1:x+1,slc,:))) );
                            [T1, T1star, A, B, res]=  fitT1_local2(Dat,TI');                      
                            T1Map(y,x,slc)=T1;
                            T1starMap(y,x,slc)=T1star;
                            AMap(y,x,slc)=A;
                            BMap(y,x,slc)=B;
                            resMap(y,x,slc,:)=res;
                        else
                            T1Map(y,x,slc)=0;
                            T1starMap(y,x,slc)=0;
                            AMap(y,x,slc)=0;
                            BMap(y,x,slc)=0;
                            resMap(y,x,slc,:)=0;
                        end        
            end
            waitbar(y/size(Dcm,1),h);
        end
        close(h);
        toc
    end
     
    warning on;
end

function [T1, T1star]=  fitT1_local(variables,TI_vect) % [S0 S1 S2 S3 ... ] [b-val0 b-val1 b-val2 .... ]
   
  
  T1=0;
  A=0;
  B=0;
  T1star=0;
  
  lb = [0 0 0];               % Lower bounds
  ub = [2000 3000 3000];            % Uper bounds
  val0 = [100 1500 1500];  % Init values  
  
  options = optimset('Display','off');
  %S0{1-2*Exp(-TI/T1)} 
  %T1_M1_relax = @(x,xdata)(  x(1)-x(2)*exp(-xdata/x(3)) ); % S0 * ( (1-f)*exp(-b*D) + f*exp(-b*Dstar) )
  
  T1_M1_relax = @(x,xdata)(  x(1)*(1-x(2)*exp(-xdata/x(3))) );
  
  [val,resnorm,residual,exitflag] = lsqcurvefit(T1_M1_relax, val0,TI_vect, variables, lb, ub, options); % lsqcurvefit(fun,x0,xdata,ydata,lb,ub,options)
  
 % S0=val(1);
% T1=val(2);
  A=val(1);
  B=val(2);
  T1star=val(3);
  T1=T1star*(B/A-1);

end


function [T1starEst, T1Est, bEst, aEst, res] = fitT1_local2(data, TI_vect)
    
%% Set the initial values for the search
[tmp,order] = sort(TI_vect);
try
  x0(3) = TI_vect(end);
  x0(1) = data(order(end));
  x0(2) = -2*x0(1);
catch
  x0 = extra.x0;
end


%% Do the fit
x = fminsearch( @(x)sum(abs( data-( x(1) + x(2)*exp(-TI_vect/x(3)) ) ).^2), x0,optimset('display','off'));

aEst = x(1) ;
bEst = x(2) ;
T1starEst = x(3);

T1Est=T1starEst*(bEst/aEst-1);
%% Compute the residual
modelValue = aEst + bEst*exp(-TI_vect/T1Est); 
res = 1/sqrt(length(data))*norm(1 - modelValue./data);

end