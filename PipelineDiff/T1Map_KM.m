

function[T1Map, S0Map]= T1Map_KM(Dcm,TI)

    T1Map=[];
    S0Map= [];
    TR=4000;
    
    

  
    lb = [0 0];            % Lower bounds
    ub = [4000 4000];            % Uper bounds
    val0 = [1000 1500];    % Init values  
  
  options = optimset('Display','off');

  T1_relax = @(x,xdata)abs(x(2)*( 1-2*exp(-(xdata/x(1)))+exp(-TR/x(1))));
  
  disp('T1 Fitting') 
  
  h = waitbar(0,'T1 fit');
  for y=1:1:size(Dcm,1)
      
      for x=1:1:size(Dcm,2)
          if Dcm(y,x,1)>0
            [val,resnorm,residual,exitflag] = lsqcurvefit(T1_relax, val0,TI, squeeze(Dcm(y,x,:))', lb, ub, options); % lsqcurvefit(fun,x0,xdata,ydata,lb,ub,options)
            T1Map(y,x)=val(1);
            S0Map(y,x)=val(2);
          else
            T1Map(y,x)=0;
            S0Map(y,x)=0;
          end
          waitbar(y/size(Dcm,1),h);
      end
  end
  close(h)
         
 end